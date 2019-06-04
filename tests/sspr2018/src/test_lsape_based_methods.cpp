/***************************************************************************
*                                                                          *
*   Copyright (C) 2018 by David B. Blumenthal                              *
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
 * @file tests/sspr2018/src/test_lsape_based_methods.cpp
 * @brief Tests performance of derived classes of ged::LSAPEBasedMethod.
 * @details The binary built from this file was used for the experiments in the following paper:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun:
 *   &ldquo;Ring based approximation of graph edit distance&rdquo;,
 *   S+SSPR 2018
 * To reproduce the experiments, install GEDLIB, go to the folder `<GEDLIB_ROOT>/tests/sspr2018/bin/` and execute the following commands:
 * @code{.sh}
 * $ cd bin
 * $ ./learn_ring_params
 * $ ./learn_subgraph_depths
 * $ ./learn_walks_depths
 * $ ./test_lsape_based_methods
 * @endcode
 * After having executed these commands, the results of the experiments are contained in the folder `<GEDLIB_ROOT>/tests/sspr2018/output/`.
 */

#define GXL_GEDLIB_SHARED
#include "../../../src/env/ged_env.hpp"

constexpr std::size_t centrality() {return 20;}
constexpr std::size_t no_centrality() {return 0;}
constexpr std::size_t lsape_optimal() {return 0;}
constexpr std::size_t gamma() {return 40;}
constexpr std::size_t lsape_greedy() {return 80;}
constexpr std::size_t no_led_option() {return 0;}

template <typename ScalarT>
void init_result_map(std::map<std::size_t, ScalarT> & result_map, ScalarT val) {
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::BIPARTITE)] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::BIPARTITE) + centrality()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::BRANCH)] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::BRANCH) + centrality()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::BRANCH_FAST)] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::BRANCH_FAST) + centrality()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::WALKS)] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::WALKS) + centrality()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::SUBGRAPH)] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::SUBGRAPH) + centrality()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + lsape_optimal()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + lsape_optimal() + centrality()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + gamma()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + gamma() + centrality()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + lsape_greedy()] = val;
	result_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + lsape_greedy() + centrality()] = val;
}

void init_method_name_map(std::map<std::size_t, std::string> & method_name_map) {
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::BIPARTITE)] = "bipartite";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::BIPARTITE) + centrality()] = "bipartite_pagerank";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::BRANCH)] = "branch";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::BRANCH) + centrality()] = "branch_pagerank";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::BRANCH_FAST)] = "branch_fast";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::BRANCH_FAST) + centrality()] = "branch_fast_pagerank";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::WALKS)] = "walks";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::WALKS) + centrality()] = "walks_pagerank";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::SUBGRAPH)] = "subgraph";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::SUBGRAPH) + centrality()] = "subgraph_pagerank";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + lsape_optimal()] = "ring_lsape_optimal";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + lsape_optimal() + centrality()] = "ring_lsape_optimal_pagerank";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + gamma()] = "ring_multisets";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + gamma() + centrality()] = "ring_multisets_pagerank";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + lsape_greedy()] = "ring_lsape_greedy";
	method_name_map[static_cast<std::size_t>(ged::Options::GEDMethod::RING) + lsape_greedy() + centrality()] = "ring_lsape_greedy_pagerank";
}

void normalize_result_map(std::map<std::size_t, double> & result_map, std::size_t val) {
	for (auto & kv : result_map) {
		kv.second /= static_cast<double>(val);
	}
}

void run_single_test(ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env, ged::GEDGraph::GraphID g_id, ged::GEDGraph::GraphID h_id,
		ged::Options::GEDMethod method, std::size_t centrality, std::size_t gamma, const std::string & options,
		std::map<std::size_t, double> & avg_time_in_sec, std::map<std::size_t, double> & avg_ub, double & best_ub, std::map<std::size_t, double> & ub,
		std::map<std::size_t, double> & distance_to_closest_graph, std::map<std::size_t, ged::GEDGraph::GraphID> & closest_graph) {
	env.set_method(method, options);
	env.run_method(g_id, h_id);
	std::size_t method_index{static_cast<std::size_t>(method) + centrality + gamma};
	double time{env.get_runtime(g_id, h_id)};
	ub[method_index] = env.get_upper_bound(g_id, h_id);
	avg_time_in_sec[method_index] += time;
	avg_ub[method_index] += ub.at(method_index);
	best_ub = std::min(best_ub, ub.at(method_index));
	if (ub.at(method_index) < distance_to_closest_graph.at(method_index)) {
		distance_to_closest_graph[method_index] = ub.at(method_index);
		closest_graph[method_index] = h_id;
	}
}

void run_tests_on_dataset(const std::string & dataset) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(std::string("../../../data/datasets/") + dataset + "/", std::string("../../../data/collections/") + dataset + ".xml"));
	if (dataset == "GREC") {
		env.set_edit_costs(ged::Options::EditCosts::GREC_1);
	}
	else {
		env.set_edit_costs(ged::Options::EditCosts::CHEM_1);
	}
	env.init();

	std::string output_prefix("../output/");
	output_prefix += dataset + "_";

	// Define option strings.
	std::string threads_option("--threads 5");
	std::string centrality_option("--centrality-method PAGERANK");
	std::string walks_depth_option("--load ");
	walks_depth_option += output_prefix + "walks.ini";
	std::string subgraph_depth_option("--load ");
	subgraph_depth_option += output_prefix + "subgraph.ini";
	std::string ring_option_lsape_optimal("--led-method LSAPE_OPTIMAL --load ");
	ring_option_lsape_optimal += output_prefix + "ring_LSAPE_OPTIMAL.ini";
	std::string ring_option_lsape_greedy("--led-method LSAPE_GREEDY --load ");
	ring_option_lsape_greedy += output_prefix + "ring_LSAPE_GREEDY.ini";
	std::string ring_option_gamma("--led-method GAMMA --load ");
	ring_option_gamma += output_prefix + "ring_GAMMA.ini";


	// Initialize result variables.
	std::map<std::size_t, double> avg_time_in_sec;
	init_result_map(avg_time_in_sec, 0.0);
	std::map<std::size_t, double> avg_ub;
	init_result_map(avg_ub, 0.0);
	std::map<std::size_t, double> classification_ratio;
	init_result_map(classification_ratio, 0.0);
	std::map<std::size_t, double> avg_dev_best_ub;
	init_result_map(avg_dev_best_ub, 0.0);
	std::map<std::size_t, double> distance_to_closest_graph;
	std::map<std::size_t, ged::GEDGraph::GraphID> closest_graph;
	std::map<std::size_t, double> ub;
	double best_ub;

	// Initialize outfile and progress bar.
	std::ofstream outfile(output_prefix + "results.csv");
	outfile << "method,avg_time_in_sec,avg_ub,classification_ratio,avg_dev_best_ub\n";
	std::size_t num_runs{(graph_ids.size() * graph_ids.size()) - graph_ids.size()};
	ged::ProgressBar progress_bar(num_runs);
	std::cout << "\r\tRunning the tests: " << progress_bar << std::flush;
	std::size_t num_best_ub_zero{0};
	for (auto g_id : graph_ids) {
		// Reset result variables used for computing classification ratio.
		init_result_map(distance_to_closest_graph, std::numeric_limits<double>::max());
		init_result_map(closest_graph, std::numeric_limits<ged::GEDGraph::GraphID>::max());
		for (auto h_id : graph_ids) {
			if (g_id == h_id) {
				continue;
			}

			// Reset result variables for computing average deviation from best upper bound.
			best_ub = std::numeric_limits<double>::max();
			init_result_map(ub, std::numeric_limits<double>::max());

			// Run methods and update result variables for computing average upper bound and runtime.
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::BIPARTITE, no_centrality(), no_led_option(), threads_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::BIPARTITE, centrality(), no_led_option(), threads_option + " " + centrality_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::BRANCH, no_centrality(), no_led_option(), threads_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::BRANCH, centrality(), no_led_option(), threads_option + " " + centrality_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::BRANCH_FAST, no_centrality(), no_led_option(), threads_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::BRANCH_FAST, centrality(), no_led_option(), threads_option + " " + centrality_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::WALKS, no_centrality(), no_led_option(), walks_depth_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::WALKS, centrality(), no_led_option(), centrality_option + " " + walks_depth_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::SUBGRAPH, no_centrality(), no_led_option(), threads_option + " " + subgraph_depth_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::SUBGRAPH, centrality(), no_led_option(), threads_option + " " + centrality_option + " " + subgraph_depth_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::RING, no_centrality(), lsape_optimal(), threads_option + " " + ring_option_lsape_optimal, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::RING, centrality(), lsape_optimal(), threads_option + " " + ring_option_lsape_optimal + " " + centrality_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::RING, no_centrality(), gamma(), threads_option + " " + ring_option_gamma, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::RING, centrality(), gamma(), threads_option + " " + ring_option_gamma + " " + centrality_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::RING, no_centrality(), lsape_greedy(), threads_option + " " + ring_option_lsape_greedy, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);
			run_single_test(env, g_id, h_id, ged::Options::GEDMethod::RING, centrality(), lsape_greedy(), threads_option + " " + ring_option_lsape_greedy + " " + centrality_option, avg_time_in_sec, avg_ub, best_ub, ub, distance_to_closest_graph, closest_graph);

			// Update result variables for computing average deviation from best upper bound.
			for (auto & kv : avg_dev_best_ub) {
				if (best_ub > 0) {
					kv.second += (ub.at(kv.first) - best_ub) / best_ub;
				}
				else {
					num_best_ub_zero++;
				}
			}

			// Update and output progress.
			progress_bar.increment();
			std::cout << "\r\tRunning the tests: " << progress_bar << std::flush;
		}

		// Update result variables for computing classification ratio.
		for (auto & kv : classification_ratio) {
			kv.second += env.get_graph_class(closest_graph.at(kv.first)) == env.get_graph_class(g_id) ? 1.0 : 0.0;
		}
	}
	std::cout << "\n";

	// Normalize the result variables.
	normalize_result_map(avg_ub, num_runs);
	normalize_result_map(avg_time_in_sec, num_runs);
	normalize_result_map(avg_dev_best_ub, num_runs - num_best_ub_zero);
	normalize_result_map(classification_ratio, graph_ids.size());

	// Write the results to the output file.
	std::map<std::size_t, std::string> method_names;
	init_method_name_map(method_names);
	for (auto & kv : method_names) {
		outfile << kv.second << "," << avg_time_in_sec.at(kv.first) << "," << avg_ub.at(kv.first) << "," << classification_ratio.at(kv.first) << "," << avg_dev_best_ub.at(kv.first) << "\n";
	}
	outfile.close();
}

int main(int argc, char* argv[]) {
	std::vector<std::string> datasets{"mao","pah","alkane","acyclic"};
	for (auto dataset : datasets) {
		try {
			run_tests_on_dataset(dataset);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << "\n";
		}
	}
	return 0;
}
