/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2019 by David B. Blumenthal                              *
 *                                                                          *
 *   This file is part of GEDLIB.                                           *
 *                                                                          *
 *   GEDLIB is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 3 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   GEDLIB is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                          *
 ***************************************************************************/

/*!
 * @file graph_clustering_heuristic.ipp
 * @brief ged::GraphClusteringHeuristic class definition.
 */

#ifndef MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTIC_IPP_
#define MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTIC_IPP_

namespace ged {

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
GraphClusteringHeuristic(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge):
ged_env_{ged_env},
mge_{mge},
main_method_{Options::GEDMethod::BRANCH_FAST},
main_options_(""),
refine_method_{Options::GEDMethod::IPFP},
refine_options_(""),
focal_graphs_("K-MEDIANS"),
init_type_("K-MEANS++"),
use_real_randomness_{true},
num_random_inits_{10},
seed_{0},
refine_{true},
time_limit_in_sec_{0},
epsilon_{0.0001},
max_itrs_{100},
print_to_stdout_{2},
assigned_focal_graph_ids_(),
node_maps_from_assigned_focal_graphs_(),
cluster_sums_of_distances_(),
cluster_radii_(),
sum_of_distances_{0},
runtime_(),
itrs_() {
	if (ged_env_ == nullptr) {
		throw Error("The environment pointer passed to the constructor of ged::GraphClusteringHeuristic is null.");
	}
	else if (not ged_env_->initialized()) {
		throw Error("The environment is uninitialized. Call ged::GEDEnv::init() before passing it to the constructor of ged::GraphClusteringHeuristic.");
	}
	else if (mge_ != nullptr) {
		if (ged_env_ != mge_->get_ged_env()) {
			throw Error("The environment pointers passed to the constructors of ged::GraphClusteringHeuristic and ged::MedianGraphEstimators don't coincide.");
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_options(const std::string & options) {
	set_default_options_();
	std::map<std::string, std::string> options_map;
	util::options_string_to_options_map(options, options_map);
	for (const auto & option : options_map) {
		if (option.first == "focal-graphs") {
			if (option.second == "MEDIANS") {
				focal_graphs_ = "K-MEDIANS";
			}
			else if (option.second == "MEDOIDS") {
				focal_graphs_ = "K-MEDOIDS";
			}
			else {
				throw ged::Error(std::string("Invalid argument ") + option.second + " for option focal-graphs. Usage: options = \"[--focal-graphs MEDIANS|MEDOIDS] [...]\"");
			}
		}
		else if (option.first == "init-type") {
			init_type_ = option.second;
			if (option.second != "K-MEANS++" and option.second != "CLUSTERS") {
				throw ged::Error(std::string("Invalid argument ") + option.second + " for option init-type. Usage: options = \"[--init-type K-MEANS++|CLUSTERS] [...]\"");
			}
		}
		else if (option.first == "random-inits") {
			try {
				num_random_inits_ = std::stoul(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option random-inits. Usage: options = \"[--random-inits <convertible to int greater 0>]\"");
			}
			if (num_random_inits_ <= 0) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option random-inits. Usage: options = \"[--random-inits <convertible to int greater 0>]\"");
			}
		}
		else if (option.first == "randomness") {
			if (option.second == "PSEUDO") {
				use_real_randomness_ = false;
			}
			else if (option.second == "REAL") {
				use_real_randomness_ = true;
			}
			else {
				throw Error(std::string("Invalid argument \"") + option.second  + "\" for option randomness. Usage: options = \"[--randomness REAL|PSEUDO] [...]\"");
			}
		}
		else if (option.first == "seed") {
			try {
				seed_ = std::stoul(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option seed. Usage: options = \"[--seed <convertible to int greater equal 0>] [...]");
			}
		}
		else if (option.first == "refine") {
			if (option.second == "TRUE") {
				refine_ = true;
			}
			else if (option.second == "FALSE") {
				refine_ = false;
			}
			else {
				throw Error(std::string("Invalid argument \"") + option.second  + "\" for option refine. Usage: options = \"[--refine TRUE|FALSE] [...]\"");
			}
		}
		else if (option.first == "time-limit") {
			try {
				time_limit_in_sec_ = std::stod(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option time-limit.  Usage: options = \"[--time-limit <convertible to double>] [...]");
			}
		}
		else if (option.first == "epsilon") {
			try {
				epsilon_ = std::stod(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option epsilon. Usage: options = \"[--epsilon <convertible to double greater 0>] [...]");
			}
			if (epsilon_ <= 0) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option epsilon. Usage: options = \"[--epsilon <convertible to double greater 0>] [...]");
			}
		}
		else if (option.first == "max-itrs") {
			try {
				max_itrs_ = std::stoi(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option max-itrs. Usage: options = \"[--max-itrs <convertible to int>] [...]");
			}
		}
		else if (option.first == "stdout") {
			if (option.second == "0") {
				print_to_stdout_ = 0;
			}
			else if (option.second == "1") {
				print_to_stdout_ = 1;
			}
			else if (option.second == "2") {
				print_to_stdout_ = 2;
			}
			else {
				throw Error(std::string("Invalid argument \"") + option.second  + "\" for option stdout. Usage: options = \"[--stdout 0|1|2] [...]\"");
			}
		}
		else {
			std::string valid_options("[--focal-graphs <arg>] [--init-type <arg>] ");
			valid_options += "[--random-inits <arg>] [--randomness <arg>] [--seed <arg>] [--refine <arg>] [--stdout <arg>] ";
			valid_options += "[--time-limit <arg>] [--max-itrs <arg>] [--epsilon <arg>] ";
			throw Error(std::string("Invalid option \"") + option.first + "\". Usage: options = \"" + valid_options + "\"");
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_main_method(Options::GEDMethod main_method, const std::string & main_options) {
	main_method_ = main_method;
	main_options_ = main_options;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_refine_method(Options::GEDMethod refine_method, const std::string & refine_options) {
	refine_method_ = refine_method;
	refine_options_ = refine_options;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
run(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids) {

	// Print information.
	if (print_to_stdout_ == 2) {
		std::cout << "\n===========================================================\n";
		std::cout << "Initializing graph clustering heuristic: ... " << std::flush;
	}

	// Sanity checks.
	if (graph_ids.empty()) {
		throw Error("Empty vector of graph IDs, unable to compute clusters.");
	}
	if (focal_graph_ids.empty()) {
		throw Error("Vector of focal graph IDs is empty, unable to compute clusters.");
	}
	if ((mge_ == nullptr) and (focal_graphs_ == "K-MEDIANS") and (max_itrs_ > 0)) {
		throw Error("No median graph estimator availabel, unable to compute clusters.");
	}

	// Start timer and record start time.
	auto start = std::chrono::high_resolution_clock::now();
	Timer timer(time_limit_in_sec_);

	// Reset information about iterations and number of times the median decreases and increases.
	itrs_ = std::vector<std::size_t>(num_random_inits_, 0);

	// Seed the random number generator.
	std::mt19937 urng;
	if (use_real_randomness_) {
		std::random_device rng;
		urng.seed(rng());
	}
	else {
		urng.seed(seed_);
	}

	// Select the employed GED method.
	ged_env_->set_method(main_method_, main_options_);

	// Initialize the best clustering.
	double best_sum_of_distances{std::numeric_limits<double>::infinity()};
	std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> best_focal_graphs;
	std::map<GEDGraph::GraphID, GEDGraph::GraphID> best_assigned_focal_graph_ids;
	std::map<GEDGraph::GraphID, NodeMap> best_node_maps_from_assigned_focal_graphs;
	std::map<GEDGraph::GraphID, double> best_distances_from_assigned_focal_graphs;
	std::map<GEDGraph::GraphID, double> best_cluster_sums_of_distances;

	// Print information.
	if (print_to_stdout_ == 2) {
		std::cout << "done.\n";
		std::cout << "-----------------------------------------------------------\n";
	}

	// Cluster the graphs.
	for (std::size_t random_init{0}; random_init < num_random_inits_; random_init++) {

		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			std::cout << "\n===========================================================\n";
			std::cout << "Running Lloyd's algorithm from initial solution " << random_init + 1 << " of " << num_random_inits_ << ".\n";
			std::cout << "-----------------------------------------------------------\n";
		}

		// Generate the initial focal graphs and the initial clusters.
		std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> focal_graphs;
		initialize_focal_graphs_and_clusters_(graph_ids, focal_graph_ids, urng, focal_graphs);

		// Run Lloyd's algorithm.
		bool converged{false};
		for (; not termination_criterion_met_(converged, timer, itrs_.at(random_init)); itrs_[random_init]++) {

			double old_sum_of_distances{sum_of_distances_};

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				std::cout << "\n===========================================================\n";
				std::cout << "Iteration " << itrs_.at(random_init) + 1 << " for initial solution " << random_init + 1 << " of " << num_random_inits_ << ".\n";
				std::cout << "-----------------------------------------------------------\n";
			}

			update_focal_graphs_(focal_graph_ids, focal_graphs);

			converged = not update_clusters_(graph_ids, focal_graph_ids, false);

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				std::cout << "Old local SOD: " << old_sum_of_distances << "\n";
				std::cout << "New local SOD: " << sum_of_distances_ << "\n";
				std::cout << "Best converged SOD: " << best_sum_of_distances << "\n";
				std::cout << "Modified clusters: " << not converged << "\n";
				std::cout << "===========================================================\n";
			}


		}

		// Update the best clustering.
		if (sum_of_distances_ < best_sum_of_distances - epsilon_) {
			best_focal_graphs = focal_graphs;
			best_sum_of_distances = sum_of_distances_;
			best_assigned_focal_graph_ids = assigned_focal_graph_ids_;
			best_node_maps_from_assigned_focal_graphs = node_maps_from_assigned_focal_graphs_;
			best_cluster_sums_of_distances = cluster_sums_of_distances_;
		}
	}

	// Print information.
	if (print_to_stdout_ == 2) {
		std::cout << "Saving the results: ... " << std::flush;
	}

	// Store the best encountered clustering.
	for (const auto & key_val : best_focal_graphs) {
		ged_env_->load_exchange_graph(key_val.second, key_val.first);
	}
	ged_env_->init(ged_env_->get_init_type());
	sum_of_distances_ = best_sum_of_distances;
	assigned_focal_graph_ids_ = best_assigned_focal_graph_ids;
	node_maps_from_assigned_focal_graphs_ = best_node_maps_from_assigned_focal_graphs;
	cluster_sums_of_distances_ = best_cluster_sums_of_distances;

	double converged_sod{sum_of_distances_};
	if (refine_) {
		update_clusters_(graph_ids, focal_graph_ids, true);
	}

	// Compute the cluster radii.
	compute_cluster_radii_(focal_graph_ids);

	// Record the runtime.
	auto end = std::chrono::high_resolution_clock::now();
	runtime_ = end - start;

	// Print global information.
	if (print_to_stdout_ != 0) {
		std::cout << "\n===========================================================\n";
		std::cout << "Finished graph clustering.\n";
		std::cout << "-----------------------------------------------------------\n";
		std::cout << "Converged SOD: " << converged_sod << "\n";
		if (refine_) {
			std::cout << "Refined SOD: " << sum_of_distances_ << "\n";
		}
		std::cout << "Cluster SODs:";
		for (const auto & key_val : cluster_sums_of_distances_) {
			std::cout << " " << key_val.second;
		}
		std::cout << "\n";
		std::cout << "Cluster radii:";
		for (const auto & key_val : cluster_radii_) {
			std::cout << " " << key_val.second;
		}
		std::cout << "\n";
		std::cout << "Overall runtime: " << runtime_.count() << "\n";
		std::cout << "Number of initial solutions: " << num_random_inits_ << "\n";
		std::size_t total_itr{0};
		for (auto itr : itrs_) {
			total_itr += itr;
		}
		std::cout << "Size of graph collection: " << graph_ids.size() << "\n";
		std::cout << "Number of clusters: " << focal_graph_ids.size() << "\n";
		std::cout << "Mean number of iterations: " << static_cast<double>(total_itr) / static_cast<double>(num_random_inits_) << "\n";
		std::cout << "===========================================================\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_runtime() const {
	return runtime_.count();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_clustering(std::map<GEDGraph::GraphID, std::vector<GEDGraph::GraphID>> & clustering) const {
	clustering.clear();
	for (const auto & key_val : cluster_radii_) {
		clustering.emplace(key_val.first, std::vector<GEDGraph::GraphID>());
	}
	for (const auto & key_val : assigned_focal_graph_ids_) {
		clustering.at(key_val.second).emplace_back(key_val.first);
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_sum_of_distances() const {
	return sum_of_distances_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_cluster_radius(GEDGraph::GraphID focal_graph_id) const {
	return cluster_radii_.at(focal_graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_cluster_sum_of_distances(GEDGraph::GraphID focal_graph_id) const {
	return cluster_sums_of_distances_.at(focal_graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDGraph::GraphID
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_assigned_focal_graph_id(GEDGraph::GraphID graph_id) const {
	return assigned_focal_graph_ids_.at(graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_distance_from_assigned_focal_graph(GEDGraph::GraphID graph_id) const {
	return node_maps_from_assigned_focal_graphs_.at(graph_id).induced_cost();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const NodeMap &
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_node_map_from_assigned_focal_graph(GEDGraph::GraphID graph_id) const {
	return node_maps_from_assigned_focal_graphs_.at(graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const std::vector<std::size_t> &
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_itrs() const {
	return itrs_;
}

template<>
void
GraphClusteringHeuristic<GXLNodeID, GXLLabel, GXLLabel>::
save(const std::string & collection_file_name, const std::string & focal_graph_dir) const {
	std::vector<std::string> gxl_file_names;
	std::vector<std::string> graph_classes;
	for (const auto & key_val : cluster_radii_) {
		GEDGraph::GraphID focal_graph_id{key_val.first};
		std::string focal_graph_name{ged_env_->get_graph_name(focal_graph_id)};
		gxl_file_names.emplace_back(focal_graph_name);
		graph_classes.emplace_back(ged_env_->get_graph_class(focal_graph_id));
		ged_env_->save_as_gxl_graph(focal_graph_id, focal_graph_dir + "/" + focal_graph_name);
	}
	ged_env_->save_graph_collection(collection_file_name, gxl_file_names, graph_classes);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_adjusted_rand_index(const std::vector<std::vector<GEDGraph::GraphID>> & ground_truth_clustering) const {

	// Ensure that the same graph IDs are contained in both clusters and transform ground truth clustering to the right format.
	std::vector<GEDGraph::GraphID> graph_ids;
	std::map<GEDGraph::GraphID, bool> contained_in_ground_truth_clustering;
	for (const auto & key_val : assigned_focal_graph_ids_) {
		contained_in_ground_truth_clustering.emplace(key_val.first, false);
		graph_ids.emplace_back(key_val.first);
	}
	std::map<GEDGraph::GraphID, std::size_t> ground_truth_cluster_ids;
	std::size_t ground_truth_cluster_id{0};
	for (const auto & ground_truth_cluster : ground_truth_clustering) {
		for (GEDGraph::GraphID graph_id : ground_truth_cluster) {
			if (contained_in_ground_truth_clustering.find(graph_id) == contained_in_ground_truth_clustering.end()) {
				throw Error("The graph with ID " + std::to_string(graph_id) + " is contained in the ground truth clustering but not in the clustered graph collection.");
			}
			contained_in_ground_truth_clustering[graph_id] = true;
			ground_truth_cluster_ids.emplace(graph_id, ground_truth_cluster_id);
		}
		ground_truth_cluster_id++;
	}
	for (const auto & key_val : contained_in_ground_truth_clustering) {
		if (not key_val.second) {
			throw Error("The graph with ID " + std::to_string(key_val.first) + " is contained in the clustered graph collection but not in the ground truth clustering.");
		}
	}

	// Compute counts of matched and un-matched pairs.
	double n_0_0{0};
	double n_0_1{0};
	double n_1_0{0};
	double n_1_1{0};
	for (auto g_id = graph_ids.begin(); g_id != graph_ids.end(); g_id++) {
		GEDGraph::GraphID cluster_id_g{assigned_focal_graph_ids_.at(*g_id)};
		std::size_t ground_truth_cluster_id_g{ground_truth_cluster_ids.at(*g_id)};
		for (auto h_id = g_id + 1; h_id != graph_ids.end(); h_id++) {
			GEDGraph::GraphID cluster_id_h{assigned_focal_graph_ids_.at(*h_id)};
			std::size_t ground_truth_cluster_id_h{ground_truth_cluster_ids.at(*h_id)};
			if (cluster_id_g != cluster_id_h) {
				if (ground_truth_cluster_id_g != ground_truth_cluster_id_h) {
					n_0_0 += 1;
				}
				else {
					n_0_1 += 1;
				}
			}
			else {
				if (ground_truth_cluster_id_g != ground_truth_cluster_id_h) {
					n_1_0 += 1;
				}
				else {
					n_1_1 += 1;
				}
			}
		}
	}

	// Return the adjusted Rand index.
	return (2 * (n_0_0 * n_1_1 - n_0_1 * n_1_0)) / ((n_0_0 + n_0_1) * (n_0_1 + n_1_1) + (n_0_0 + n_1_0) * (n_1_0 + n_1_1));

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_gini_coefficient() const {

	// Compute cluster sizes.
	std::map<GEDGraph::GraphID, int> cluster_sizes;
	for (const auto & key_val : cluster_radii_) {
		cluster_sizes.emplace(key_val.first, 0);
	}
	for (const auto & key_val : assigned_focal_graph_ids_) {
		cluster_sizes[key_val.second]++;
	}

	// Compute the numerator.
	int numerator{0};
	for (const auto & key_val_1 : cluster_sizes) {
		for (const auto & key_val_2 : cluster_sizes) {
			numerator += std::abs(key_val_1.second - key_val_2.second);
		}
	}

	// Return the Gini coefficient.
	return static_cast<double>(numerator) / static_cast<double>(cluster_sizes.size() * assigned_focal_graph_ids_.size());

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_silhouette_score() const {
    if (refine_) {
        ged_env_->set_method(refine_method_, refine_options_);
    }
    else if (focal_graphs_ == "K-MEDIANS") {
        ged_env_->set_method(main_method_, main_options_);
    }
    double score{0};
    for (const auto & key_val : assigned_focal_graph_ids_) {
        GEDGraph::GraphID g_id{key_val.first};
        GEDGraph::GraphID assigned_focal_graph_id{key_val.second};
        std::map<GEDGraph::GraphID, double> mean_cluster_dists;
        std::map<GEDGraph::GraphID, std::size_t> cluster_sizes;
        for (const auto & key_val_2 : cluster_radii_) {
            mean_cluster_dists.emplace(key_val_2.first, 0);
            cluster_sizes.emplace(key_val_2.first, 0);
        }
        for (const auto & key_val_2 : assigned_focal_graph_ids_) {
            GEDGraph::GraphID h_id{key_val_2.first};
            GEDGraph::GraphID focal_graph_id{key_val_2.second};
            if (g_id == h_id) {
                continue;
            }
            cluster_sizes[focal_graph_id]++;
            ged_env_->run_method(g_id, h_id);
            mean_cluster_dists[focal_graph_id] += ged_env_->get_upper_bound(g_id, h_id);
        }
        if (cluster_sizes.at(assigned_focal_graph_id) == 0) {
            continue;
        }
        double a{0};
        double b{std::numeric_limits<double>::infinity()};
        for (auto & key_val_2 : mean_cluster_dists) {
            GEDGraph::GraphID focal_graph_id{key_val_2.first};
            key_val_2.second /= static_cast<double>(cluster_sizes.at(focal_graph_id));
            if (focal_graph_id == assigned_focal_graph_id) {
                a = key_val_2.second;
            }
            else if (key_val_2.second < b) {
                b = key_val_2.second;
            }
        }
        score += (b - a) / std::max(b, a);
    }
    return score / static_cast<double>(assigned_focal_graph_ids_.size());
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> *
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_ged_env() {
	return ged_env_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_default_options_() {
	focal_graphs_ = "K-MEDIANS";
	init_type_ = "K-MEANS++";
	num_random_inits_ = 10;
	use_real_randomness_ = true;
	seed_ = 0;
	time_limit_in_sec_ = 0;
	epsilon_ = 0.0001;
	max_itrs_ = 100;
	print_to_stdout_ = 2;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
initialize_focal_graphs_and_clusters_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
		std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs) {

	// Generate the focal graphs.
	focal_graphs.clear();
	if (init_type_ == "K-MEANS++") {
		compute_initial_focal_graphs_k_means_plus_plus_(graph_ids, focal_graph_ids, urng, focal_graphs);
	}
	else {
		compute_initial_focal_graphs_cluster_sampling_(graph_ids, focal_graph_ids, urng, focal_graphs);
	}

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "Initializing focal graphs: ... " << std::flush;
	}

	// Initialize the clusters.
	assigned_focal_graph_ids_.clear();
	node_maps_from_assigned_focal_graphs_.clear();
	cluster_sums_of_distances_.clear();
	NodeMap empty_node_map(0, 0);
	for (GEDGraph::GraphID graph_id : graph_ids) {
		assigned_focal_graph_ids_.emplace(graph_id, undefined());
		node_maps_from_assigned_focal_graphs_.emplace(graph_id, empty_node_map);
	}
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
		cluster_sums_of_distances_.emplace(focal_graph_id, std::numeric_limits<double>::infinity());
		cluster_radii_.emplace(focal_graph_id, std::numeric_limits<double>::infinity());
	}
	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "done.\n";
	}

	update_clusters_(graph_ids, focal_graph_ids, false);

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_initial_focal_graphs_k_means_plus_plus_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
		std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs) {

	// Print information about current iteration.
	ged::ProgressBar progress(focal_graph_ids.size());
	if (print_to_stdout_ == 2) {
		std::cout << "\rk-means++ initialization: " << progress << std::flush;
	}

	// IDs of the graphs that will be selected as focal graphs.
	std::vector<GEDGraph::GraphID> original_focal_graph_ids;

	// Flags that indicate for all graphs whether they have been selected as focal graphs.
	std::vector<bool> already_selected(graph_ids.size(), false);

	// Select the first focal graph.
	std::uniform_int_distribution<std::size_t> unif_dist(0, graph_ids.size() - 1);
	GEDGraph::GraphID selected_graph_id_pos{unif_dist(urng)};
	original_focal_graph_ids.emplace_back(graph_ids.at(selected_graph_id_pos));
	focal_graphs.emplace(focal_graph_ids.at(focal_graphs.size()), ged_env_->get_graph(original_focal_graph_ids.back()));
	already_selected[selected_graph_id_pos] = true;

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		progress.increment();
		std::cout << "\rk-means++ initialization: " << progress << std::flush;
	}

	// Select all remaining focal graphs.
	while (focal_graphs.size()  < focal_graph_ids.size()) {

		// Compute the weights that yield the probability distribution.
		std::vector<double> weights(graph_ids.size(), std::numeric_limits<double>::infinity());
		for (std::size_t graph_id_pos{0}; graph_id_pos < graph_ids.size(); graph_id_pos++) {
			if (already_selected.at(graph_id_pos)) {
				weights[graph_id_pos] = 0;
				continue;
			}
			for (const auto & focal_graph_id : original_focal_graph_ids) {
				ged_env_->run_method(focal_graph_id, graph_ids.at(graph_id_pos));
				weights[graph_id_pos] = std::min(weights.at(graph_id_pos), ged_env_->get_upper_bound(focal_graph_id, graph_ids.at(graph_id_pos)));
			}
		}

		// Draw the next focal graph from the obtained distribution.
		std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
		selected_graph_id_pos = dist(urng);
		original_focal_graph_ids.emplace_back(graph_ids.at(selected_graph_id_pos));
		focal_graphs.emplace(focal_graph_ids.at(focal_graphs.size()), ged_env_->get_graph(original_focal_graph_ids.back()));
		already_selected[selected_graph_id_pos] = true;

		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			progress.increment();
			std::cout << "\rk-means++ initialization: " << progress << std::flush;
		}
	}

	// Load the focal graphs into the environment.
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
		ged_env_->load_exchange_graph(focal_graphs.at(focal_graph_id), focal_graph_id);
	}
	ged_env_->init(ged_env_->get_init_type());

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n";
	}

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_initial_focal_graphs_cluster_sampling_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
		std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs) {

	// Print information about current iteration.
	ged::ProgressBar progress(focal_graph_ids.size());
	if (print_to_stdout_ == 2) {
		std::cout << "\rcluster initialization: " << progress << std::flush;
	}

	// Determine the cluster size.
	std::size_t cluster_size{graph_ids.size() / focal_graph_ids.size()};

	// Shuffle the graph IDs.
	std::vector<GEDGraph::GraphID> shuffled_graph_ids(graph_ids);
	std::shuffle(shuffled_graph_ids.begin(), shuffled_graph_ids.end(), urng);

	// Initialize the clustering.
	std::vector<std::vector<GEDGraph::GraphID>> clustering;

	// Populate the clustering.
	std::size_t pos{0};
	while (clustering.size() < focal_graph_ids.size() - 1) {
		clustering.emplace_back();
		for (; pos < clustering.size() * cluster_size; pos++) {
			clustering.back().emplace_back(shuffled_graph_ids.at(pos));
		}
	}
	clustering.emplace_back();
	for (; pos < shuffled_graph_ids.size(); pos++) {
		clustering.back().emplace_back(shuffled_graph_ids.at(pos));
	}

	// Compute the initial focal graphs as medians or medoids of the initial clusters.
	for (std::size_t cluster_id{0}; cluster_id < clustering.size(); cluster_id++) {
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> focal_graph;
		std::map<GEDGraph::GraphID, NodeMap> new_node_maps_from_focal_graph;
		if (focal_graphs_ == "K-MEDIANS") {
			compute_median_(clustering.at(cluster_id), focal_graph_ids.at(cluster_id), focal_graph, new_node_maps_from_focal_graph);
		}
		else {
			compute_medoid_(clustering.at(cluster_id), focal_graph_ids.at(cluster_id), focal_graph, new_node_maps_from_focal_graph);
		}
		focal_graphs.emplace(focal_graph_ids.at(cluster_id), focal_graph);

		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			progress.increment();
			std::cout << "\rcluster initialization: " << progress << std::flush;
		}
	}
	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
termination_criterion_met_(bool converged, const Timer & timer, std::size_t itr) const {
	return converged or timer.expired() or (max_itrs_ >= 0 ? itr >= max_itrs_ : false);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_clusters_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids, bool refine) {

	// Print information about current iteration.
	ProgressBar progress(graph_ids.size() * focal_graph_ids.size());
	if (print_to_stdout_ == 2) {
		if (refine) {
			std::cout << "\rRefining clusters: " << progress << std::flush;
		}
		else {
			std::cout << "\rUpdating clusters: " << progress << std::flush;
		}
	}

	// Re-select the employed GED method if used during refinement or if K-MEDIANS is used (the median graph estimator might have used another method).
	if (refine) {
		ged_env_->set_method(refine_method_, refine_options_);
	}
	else if (focal_graphs_ == "K-MEDIANS") {
		ged_env_->set_method(main_method_, main_options_);
	}

	// Determine the closest focal graph for each graph.
	bool clusters_modified{false};
	for (GEDGraph::GraphID graph_id : graph_ids) {
		GEDGraph::GraphID old_assigned_focal_graph_id{assigned_focal_graph_ids_.at(graph_id)};
		for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {

			// Update the assigned focal graph.
			ged_env_->run_method(focal_graph_id, graph_id);
			if (ged_env_->get_upper_bound(focal_graph_id, graph_id) < node_maps_from_assigned_focal_graphs_.at(graph_id).induced_cost() - epsilon_) {
				assigned_focal_graph_ids_[graph_id] = focal_graph_id;
				node_maps_from_assigned_focal_graphs_.at(graph_id) = ged_env_->get_node_map(focal_graph_id, graph_id);
			}

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				progress.increment();
				if (refine) {
					std::cout << "\rRefining clusters: " << progress << std::flush;
				}
				else {
					std::cout << "\rUpdating clusters: " << progress << std::flush;
				}
			}
		}
		if (assigned_focal_graph_ids_.at(graph_id) != old_assigned_focal_graph_id) {
			clusters_modified = true;
		}
	}

	// Update the sums of distances for the clusters and the overall sum of distances.
	if (clusters_modified) {
		sum_of_distances_ = 0;
		for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
			cluster_sums_of_distances_.at(focal_graph_id) = 0;
		}
		for (const auto & key_val : assigned_focal_graph_ids_) {
			double new_cost{node_maps_from_assigned_focal_graphs_.at(key_val.first).induced_cost()};
			cluster_sums_of_distances_.at(key_val.second) += new_cost;
			sum_of_distances_ += new_cost;
		}
	}

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n";
	}


	// Return true if the clusters were modified.
	return clusters_modified;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_focal_graphs_(const std::vector<GEDGraph::GraphID> & focal_graph_ids, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs) {

	// Print information about current iteration.
	ProgressBar progress(focal_graph_ids.size());
	if (print_to_stdout_ == 2) {
		std::cout << "\rUpdating focal graphs: " << progress << std::flush;
	}

	// Collect the clusters.
	std::map<GEDGraph::GraphID, std::vector<GEDGraph::GraphID>> clustering;
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
		clustering.emplace(focal_graph_id, std::vector<GEDGraph::GraphID>());
	}
	for (const auto & key_val : assigned_focal_graph_ids_) {
		clustering.at(key_val.second).emplace_back(key_val.first);
	}

	// Update the focal graphs, the sums of distances for the clusters, and the node maps.
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {

		// Continue with the next cluster if the current cluster is empty.
		if (clustering.at(focal_graph_id).empty()) {
			continue;
			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				progress.increment();
				std::cout << "\rUpdating focal graphs: " << progress << std::flush;
			}
		}

		// Compute new focal graph for the current cluster.
		std::map<GEDGraph::GraphID, NodeMap> new_node_maps_from_focal_graph;
		double new_cluster_sum_of_distances{0};
		if (focal_graphs_ == "K-MEDIANS") {
			new_cluster_sum_of_distances = compute_median_(clustering.at(focal_graph_id), focal_graph_id, focal_graphs.at(focal_graph_id), new_node_maps_from_focal_graph);
		}
		else {
			new_cluster_sum_of_distances = compute_medoid_(clustering.at(focal_graph_id), focal_graph_id, focal_graphs.at(focal_graph_id), new_node_maps_from_focal_graph);
		}
		if (new_cluster_sum_of_distances < cluster_sums_of_distances_.at(focal_graph_id) - epsilon_) {
			cluster_sums_of_distances_.at(focal_graph_id) = new_cluster_sum_of_distances;
			for (GEDGraph::GraphID graph_id : clustering.at(focal_graph_id)) {
				node_maps_from_assigned_focal_graphs_.at(graph_id) = new_node_maps_from_focal_graph.at(graph_id);
			}
		}

		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			progress.increment();
			std::cout << "\rUpdating focal graphs: " << progress << std::flush;
		}
	}

	// Update the overall sum of distances.
	sum_of_distances_ = 0;
	for (const auto & key_val : cluster_sums_of_distances_) {
		sum_of_distances_ += key_val.second;
	}

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_median_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID median_id,
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median, std::map<GEDGraph::GraphID, NodeMap> & node_maps_from_median) const {

	// Determine the median and load it into the environment.
	mge_->run(graph_ids, median_id);

	// Update the median's ExchangeGraph representation.
	median = ged_env_->get_graph(median_id);

	// Update the node maps within the cluster.
	node_maps_from_median.clear();
	for (GEDGraph::GraphID graph_id : graph_ids) {
		node_maps_from_median.emplace(graph_id, mge_->get_node_map_from_median(graph_id));
	}

	// Return the sum of distances.
	return mge_->get_sum_of_distances();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_medoid_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID medoid_id,
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & medoid, std::map<GEDGraph::GraphID, NodeMap> & node_maps_from_medoid) const {

	// Determine the medoid.
	GEDGraph::GraphID original_medoid_id{undefined()};
	double best_sum_of_distances{std::numeric_limits<double>::infinity()};
	for (auto g_id : graph_ids) {
		double sum_of_distances{0};
		for (auto h_id : graph_ids) {
			ged_env_->run_method(g_id, h_id);
			sum_of_distances += ged_env_->get_upper_bound(g_id, h_id);
		}
		if (sum_of_distances < best_sum_of_distances) {
			best_sum_of_distances = sum_of_distances;
			original_medoid_id = g_id;
		}
	}

	// Update the medoid's exchange graph representation.
	medoid = ged_env_->get_graph(original_medoid_id);

	// Update the node maps within the cluster.
	node_maps_from_medoid.clear();
	for (auto graph_id : graph_ids) {
		node_maps_from_medoid.emplace(graph_id, ged_env_->get_node_map(original_medoid_id, graph_id));
	}

	// Load the medoid into the environment.
	ged_env_->load_exchange_graph(medoid, medoid_id);
	ged_env_->init(ged_env_->get_init_type());

	// Return the sum of distances.
	return best_sum_of_distances;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_cluster_radii_(const std::vector<GEDGraph::GraphID> & focal_graph_ids) {

	// Print information.
	ProgressBar progress(assigned_focal_graph_ids_.size());
	if (print_to_stdout_ == 2) {
		std::cout << "\rComputing cluster radii: " << progress << std::flush;
	}

	// Initialize the cluster radii.
	cluster_radii_.clear();
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
		cluster_radii_.emplace(focal_graph_id, 0);
	}

	// Compute the radii.
	for (const auto & key_val : assigned_focal_graph_ids_) {
		double old_radius{cluster_radii_.at(key_val.second)};
		double current_cost{node_maps_from_assigned_focal_graphs_.at(key_val.first).induced_cost()};
		cluster_radii_.at(key_val.second) = std::max(old_radius, current_cost);
		// Print information.
		if (print_to_stdout_ == 2) {
			progress.increment();
			std::cout << "\rComputing cluster radii: " << progress << std::flush;
		}
	}

	// Print information.
	if (print_to_stdout_ == 2) {
		std::cout << "\n";
	}
}

}


#endif /* MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTIC_IPP_ */
