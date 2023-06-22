//
// Created by neon on 23-6-3.
//

#include <thrust/copy.h>
#include "../env/matrix.hpp"
#include "../env/ged_graph.hpp"

namespace ged {
    __device__ double *ptr{};
    __device__ std::size_t **device_sorted_edge_labels_g;
    __device__ std::size_t **device_sorted_edge_labels_h;
    __device__ int *sizes_of_array_g;
    __device__ int *sizes_of_array_h;
    __device__ double *dummy_row{}, *dummy_col{};
    __device__ int *g_deg_data{}, *h_deg_data{};
    __device__ double **relabeling_costs;
    __device__ unsigned **intersection_costs;
    __device__ double **results{};

    __device__ cudaStream_t insertion, deletion, relabeling, intersection;

    static std::size_t omp_num_threads;

    extern "C"
    __host__ void
    prepare_cuda_env_(const ged::DMatrix &mat,
                      const std::map<GEDGraph::NodeID, std::vector<LabelID>> &eg,
                      const std::map<GEDGraph::NodeID, std::vector<LabelID>> &eh,
                      vector<int> &g_deg,
                      vector<int> &h_deg,
                      std::size_t num_threads
    ) {
        omp_num_threads = num_threads;

        // 有问题 没有初始化?
        cudaMalloc((void **) &ptr, (mat.num_rows() * mat.num_rows()) * sizeof(size_t));
        cudaMalloc((void **) &dummy_row, sizeof(double) * mat.num_rows());
        cudaMalloc((void **) &dummy_col, sizeof(double) * mat.num_cols());
        thrust::copy(mat.data(), mat.data() + mat.num_rows() * mat.num_rows(), ptr);
        for (int i = 0; i < mat.num_rows(); i++) {
            dummy_row[i] = mat(0, i);
        }
        for (int i = 0; i < mat.num_cols(); i++) {
            dummy_col[i] = mat(i, 0);
        }
        // Tricky 有空修
        // 可以改异步io
        cudaMalloc(device_sorted_edge_labels_g, eg.size() * sizeof(std::size_t));
        cudaMalloc((void **) &sizes_of_array_g, eg.size() * sizeof(int));
        for (int i = 0; i < eg.size(); i++) {
            auto &p = device_sorted_edge_labels_g[i];
            //cudaMalloc((void**)&p,sizeof(std::size_t) * eg.size());
            auto item = eg.find(i)->second;
            sizes_of_array_g[i] = item.size();
            thrust::copy(item.data(), item.data() + item.size(), p);

        }
        cudaMalloc(device_sorted_edge_labels_h, eg.size() * sizeof(std::size_t));
        cudaMalloc((void **) &sizes_of_array_h, eg.size() * sizeof(int));
        for (int i = 0; i < eh.size(); i++) {
            auto &p = device_sorted_edge_labels_h[i];
            //cudaMalloc((void**)&p,sizeof(std::size_t) * eg.size());
            auto item = eh.find(i)->second;
            sizes_of_array_h[i] = item.size();
            thrust::copy(item.data(), item.data() + item.size(), p);

        }
        cudaMalloc((void **) &g_deg_data, sizeof(int) * g_deg.size());
        cudaMemcpy(g_deg_data, g_deg.data(), g_deg.size(), cudaMemcpyHostToDevice);


        cudaMalloc((void **) &h_deg_data, sizeof(int) * h_deg.size());
        cudaMemcpy(h_deg_data, h_deg.data(), h_deg.size(), cudaMemcpyHostToDevice);
    }

    __global__ void
    compute_deletion_cost(ged::GEDGraph::NodeID i,ged::GEDGraph::NodeID k,int sub) {
        extern __shared__ double reduction_helper[];

        auto tid = threadIdx.x;
        auto idx = device_sorted_edge_labels_h[k][tid];

        reduction_helper[tid] = dummy_row[idx];

        for (int s = blockDim.x / 2; s > 0; s >>= 1) {
            if (tid < s) {
                reduction_helper[tid + s] = thrust::min(reduction_helper[tid], reduction_helper[tid + s]);
            }
            __syncthreads();
        }

        if (tid == 0) {
            results[i][k] = reduction_helper[0] * sub * 0.5;
        }
    }

    __global__ void
    compute_insertion_cost(ged::GEDGraph::NodeID i,ged::GEDGraph::NodeID k,int sub) {
        extern __shared__ double reduction_helper[];

        auto tid = threadIdx.x;
        auto idx = device_sorted_edge_labels_g[i][tid];

        reduction_helper[tid] = dummy_col[idx];

        for (auto s = blockDim.x / 2; s > 0; s >>= 1) {
            if (tid < s) {
                reduction_helper[tid + s] = thrust::min(reduction_helper[tid], reduction_helper[tid + s]);
            }
            __syncthreads();
        }

        if (tid == 0) {
            results[i][k] = reduction_helper[0] * sub * 0.5;
        }
    }

    __global__ void
    compute_relabelling_cost(ged::GEDGraph::NodeID i, ged::GEDGraph::NodeID k) {
        extern __shared__ double reduction_helper[];

        // todo need recheck
        auto tid = threadIdx.x;
        auto bid = blockIdx.x;
        auto col_size = gridDim.x;
        auto pos = bid * blockDim.x + tid;
        double relabel_cost = 1e300;

        // Step1. filter out the nodes that need relabeled
        auto row = device_sorted_edge_labels_g[i][bid];
        auto col = device_sorted_edge_labels_h[k][tid];
        if (row != col)
            relabel_cost = ptr[row * col_size + col];

        // Step2. Reduce
        reduction_helper[pos] = relabel_cost;
        for (auto s = (gridDim.x * blockDim.x / 2); s > 0; s >>= 1) {
            if (pos < s) {
                reduction_helper[pos] = thrust::min(reduction_helper[pos], reduction_helper[s + pos]);
            }
            __syncthreads();
        }
        if (pos == 0) {
            relabeling_costs[i][k] = reduction_helper[0];
        }
    }


    __global__ void
    compute_multiset_intersection_size(ged::GEDGraph::NodeID i, ged::GEDGraph::NodeID k) {
        // Worst case: O(n)
        unsigned __shared__ offset_g, offset_h;
        unsigned size{};
        auto tid = threadIdx.x;
        auto sg = sizes_of_array_g[i];
        auto sh = sizes_of_array_h[k];
        offset_g = offset_h  = 0; // <- is this necessary ?

        while (true) {
            auto lg = device_sorted_edge_labels_g[i][tid + offset_g];
            auto lh = device_sorted_edge_labels_h[k][tid + offset_h];
            if (lg == lh) {
                atomicAdd(&size, 1);
                break;
            } else {
                if (tid == 0) {
                    if (lg > lh) {
                        offset_h++;
                    } else {
                        offset_g++;
                    }
                }
            }
            if (tid + offset_g >= sg || tid + offset_h >= sh)
                break;
            __syncthreads();
        }
        intersection_costs[i][k] = size;
    }

    __global__ void
    compute_substitution_cost_with_cuda(
    ) {
        // Collect node substitution cost.
        double cost{};

        auto i = blockIdx.x;
        auto k = threadIdx.x;
        auto g_deg = g_deg_data[i];
        auto h_deg = h_deg_data[k];
        /*auto g_num_nodes = gridDim.x;
        auto h_num_nodes = blockDim.x;
        auto sg = sizes_of_array_g[i];
        auto sh = sizes_of_array_h[k];*/

        double min_relabling_cost = relabeling_costs[i][k];
        unsigned intersection_size = intersection_costs[i][k];

        // Write Back
        if (thrust::min(h_deg, g_deg) - intersection_size > 0) {
            cost += static_cast<double>(thrust::min(g_deg, h_deg) - intersection_size) * min_relabling_cost * 0.5;
        }
        results[i][k] += cost;
    }

    extern "C"
    __host__  double *
    launch_kernel(int g_num_nodes, int h_num_nodes,vector<int>& g_degs,vector<int>& h_degs) {
        double **ret{};
        cudaMalloc(results,sizeof(double)* g_num_nodes * h_num_nodes);
        cudaMalloc((void**) &relabeling_costs,sizeof(double) * g_num_nodes * h_num_nodes);
        cudaMalloc(intersection_costs,sizeof(unsigned) * g_num_nodes * h_num_nodes);

        cudaStreamCreate(&insertion);
        cudaStreamCreate(&deletion);
        cudaStreamCreate(&relabeling);
        cudaStreamCreate(&intersection);

#ifdef _OPENMP
        omp_set_num_threads(omp_num_threads - 1);
#pragma omp parallel for if(omp_num_threads > 1)
#endif
        for(auto i = 0;i<g_num_nodes;i++) {
            for (auto  k = 0;k<h_num_nodes;k++) {
                auto g_deg = g_degs[i];
                auto h_deg = h_degs[k];

                // 1. 核函数是否接受__device__参数
                // 2. malloc 2层指针是否可以按照矩阵形式访问

                // Compute insertion cost.
                if (g_deg < h_deg) {
                    compute_deletion_cost
                    <<<1, sizes_of_array_h[k], sizes_of_array_h[k] * sizeof(double), insertion>>>
                            (i,k,h_deg - g_deg);
                }

                // Compute deletion cost.
                if (g_deg > h_deg) {
                    compute_insertion_cost
                    <<<1, sizes_of_array_g[i], sizes_of_array_g[i] * sizeof(double), deletion>>>
                            (i,k,g_deg - h_deg);
                }

                // Compute relabeling cost.
                compute_relabelling_cost
                <<<g_num_nodes, h_num_nodes, g_num_nodes * h_num_nodes * sizeof(double), relabeling>>>
                        (i, k);

                // Compute multiset intersection size.
                compute_multiset_intersection_size<<<1, thrust::min(sizes_of_array_g[i], sizes_of_array_h[k])>>>(i, k);
            }
        }
        cudaDeviceSynchronize();

        compute_substitution_cost_with_cuda<<<g_num_nodes,h_num_nodes>>>();
        cudaStreamDestroy(insertion);
        cudaStreamDestroy(deletion);
        cudaStreamDestroy(relabeling);
        cudaStreamDestroy(intersection);

        cudaMallocHost(ret,g_num_nodes * h_num_nodes * sizeof(double));
        cudaMemcpy(ret,results,g_num_nodes * h_num_nodes * sizeof(double),cudaMemcpyDeviceToHost);
        return nullptr;
    }
}