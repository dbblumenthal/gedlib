//
// Created by David Blumenthal on 06.10.20.
//

#define GXL_GEDLIB_SHARED
#include "util.hpp"

void test_on_dataset(const std::string & dataset, std::size_t num_pairs, std::size_t num_permutations) {

    // Initialize environment.
    std::cout << "\n=== " << dataset << " ===\n";
    std::cout << "\tInitializing the environment ...\n";
    ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
    util::setup_environment(dataset, env);
    std::vector<std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID>> all_pairs;
    for (ged::GEDGraph::GraphID g_id{0}; g_id < env.num_graphs() - 1; g_id++) {
        for (ged::GEDGraph::GraphID h_id{g_id + 1}; h_id < env.num_graphs(); h_id++) {
            all_pairs.emplace_back(g_id, h_id);
        }
    }
    auto rng = std::default_random_engine {};
    std::shuffle(all_pairs.begin(), all_pairs.end(), rng);
    std::size_t real_num_pairs{std::min(num_pairs, all_pairs.size())};
    std::vector<std::size_t> nums_solutions{1, 5, 10, 15, 20};
    std::map<std::size_t, ged::DMatrix> runtimes_baseline;
    std::map<std::size_t, ged::DMatrix> runtimes_dissimilar;
    std::map<std::size_t, ged::DMatrix> ubs_baseline;
    std::map<std::size_t, ged::DMatrix> ubs_dissimilar;
    for (auto num_solutions : nums_solutions) {
        runtimes_baseline[num_solutions] = ged::DMatrix(real_num_pairs, num_permutations);
        runtimes_dissimilar[num_solutions] = ged::DMatrix(real_num_pairs, num_permutations);
        ubs_baseline[num_solutions] = ged::DMatrix(real_num_pairs, num_permutations);
        ubs_dissimilar[num_solutions] = ged::DMatrix(real_num_pairs, num_permutations);
    }
    ged::ProgressBar progress(num_permutations * nums_solutions.size());
    std::cout << "\rRunning tests: " << progress << std::flush;
    for (std::size_t perm_id{0}; perm_id < num_permutations; perm_id++) {
        ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
        util::setup_environment(dataset, env);
        for (std::size_t num_solutions : nums_solutions) {
            std::string options("--max-num-solutions " + std::to_string(num_solutions) + " --enumeration-method ");
            for (std::size_t pair_id{0}; pair_id < real_num_pairs; pair_id++) {
                ged::GEDGraph::GraphID g_id{all_pairs.at(pair_id).first};
                ged::GEDGraph::GraphID h_id{all_pairs.at(pair_id).second};
                env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, options + "BASELINE");
                env.run_method(g_id, h_id, true);
                runtimes_baseline.at(num_solutions)(pair_id, perm_id) = env.get_runtime(g_id, h_id);
                ubs_baseline.at(num_solutions)(pair_id, perm_id) = env.get_upper_bound(g_id, h_id);
                env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, options + "DISSIMILAR");
                env.run_method(g_id, h_id, true);
                runtimes_dissimilar.at(num_solutions)(pair_id, perm_id) = env.get_runtime(g_id, h_id);
                ubs_dissimilar.at(num_solutions)(pair_id, perm_id) = env.get_upper_bound(g_id, h_id);
            }
            progress.increment();
            std::cout << "\rRunning tests: " << progress << std::flush;
        }
    }
    std::vector<double> mean_mean_runtimes_baseline;
    std::vector<double> mean_mean_runtimes_dissimilar;
    std::vector<double> mean_mean_ubs_baseline;
    std::vector<double> mean_mean_ubs_dissimilar;
    std::vector<double> improvement_mean_mean_ubs_baseline;
    std::vector<double> improvement_mean_mean_ubs_dissimilar;
    std::vector<double> mean_std_ubs_baseline;
    std::vector<double> mean_std_ubs_dissimilar;
    std::vector<double> improvement_mean_std_ubs_baseline;
    std::vector<double> improvement_mean_std_ubs_dissimilar;
    std::cout << "\naggregating results\n";
    for (auto num_solutions : nums_solutions) {
        mean_mean_runtimes_baseline.emplace_back(runtimes_baseline.at(num_solutions).matrix().rowwise().mean().mean());
        auto mean_ubs_baseline = ubs_baseline.at(num_solutions).matrix().rowwise().mean();
        mean_mean_ubs_baseline.emplace_back(mean_ubs_baseline.mean());
        auto squared_error_ubs_baseline = (ubs_baseline.at(num_solutions).matrix().colwise() - mean_ubs_baseline).array().square();
        mean_std_ubs_baseline.emplace_back((squared_error_ubs_baseline.rowwise().sum()/(num_permutations - 1)).sqrt().mean());
        mean_mean_runtimes_dissimilar.emplace_back(runtimes_dissimilar.at(num_solutions).matrix().rowwise().mean().mean());
        auto mean_ubs_dissimilar = ubs_dissimilar.at(num_solutions).matrix().rowwise().mean();
        mean_mean_ubs_dissimilar.emplace_back(mean_ubs_dissimilar.mean());
        auto squared_error_ubs_dissimilar = (ubs_dissimilar.at(num_solutions).matrix().colwise() - mean_ubs_dissimilar).array().square();
        mean_std_ubs_dissimilar.emplace_back((squared_error_ubs_dissimilar.rowwise().sum()/(num_permutations - 1)).sqrt().mean());
    }
    std::cout << "writing output\n";
    std::ofstream outfile("../output/" + dataset + ".csv");
    outfile << "num_solutions,mean_mean_runtime_baseline,mean_mean_runtime_dissimilar,";
    outfile << "mean_mean_ub_baseline,mean_mean_ub_dissimilar,mean_std_ub_baseline,mean_std_ub_dissimilar,";
    outfile << "improvement_mean_mean_ub_baseline,improvement_mean_mean_ub_dissimilar,improvement_mean_std_ub_baseline,improvement_mean_std_ub_dissimilar\n";
    for (std::size_t sol_id{0}; sol_id < nums_solutions.size(); sol_id++) {
        outfile << nums_solutions.at(sol_id) << ",";
        outfile << mean_mean_runtimes_baseline.at(sol_id) << ",";
        outfile << mean_mean_runtimes_dissimilar.at(sol_id) << ",";
        outfile << mean_mean_ubs_baseline.at(sol_id) << ",";
        outfile << mean_mean_ubs_dissimilar.at(sol_id) << ",";
        outfile << mean_std_ubs_baseline.at(sol_id) << ",";
        outfile << mean_std_ubs_dissimilar.at(sol_id) << ",";
        outfile << 100 * (mean_mean_ubs_baseline.at(0) - mean_mean_ubs_baseline.at(sol_id)) / mean_mean_ubs_baseline.at(0) << ",";
        outfile << 100 * (mean_mean_ubs_dissimilar.at(0) - mean_mean_ubs_dissimilar.at(sol_id)) / mean_mean_ubs_dissimilar.at(0) << ",";
        outfile << 100 * (mean_std_ubs_baseline.at(0) - mean_std_ubs_baseline.at(sol_id)) / mean_std_ubs_baseline.at(0) << ",";
        outfile << 100 * (mean_std_ubs_dissimilar.at(0) - mean_std_ubs_dissimilar.at(sol_id)) / mean_std_ubs_dissimilar.at(0) << "\n";
    }
    outfile.close();
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "usage: " << argv[0] << " <dataset> <num_permutations> <num_graph_pairs>" << std::endl;
        return 1;
    }
    std::string dataset(argv[1]);
    std::size_t num_permutations(atoi(argv[2]));
    std::size_t num_pairs(atoi(argv[3]));
    try {
        test_on_dataset(dataset, num_pairs, num_permutations);
        return 0;
    }
    catch (...) {
        return 1;
    }
}

