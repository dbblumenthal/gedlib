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
 * @file classification_tests.cpp
 * @brief Tests the GraphClusteringHeuristic when used for data reduction and classification.
 */

#define GXL_GEDLIB_SHARED

#include "../src/graph_bst.hpp"

std::unordered_set<std::string> irrelevant_node_attributes(const std::string & dataset) {
	std::unordered_set<std::string> irrelevant_attributes;
	if (dataset == "AIDS") {
		irrelevant_attributes.insert({"x", "y", "symbol", "charge"});
	}
	return irrelevant_attributes;
}

bool constant_node_costs(const std::string & dataset) {
	if (dataset == "Letter") {
		return false;
	}
	else if (dataset != "Mutagenicity" and dataset != "AIDS") {
		throw ged::Error("Invalid dataset " + dataset + ". Usage: ./median_tests <AIDS|Mutagenicity|Letter>");
	}
	return true;
}

ged::Options::EditCosts edit_costs(const std::string & dataset) {
	if (dataset == "Letter") {
		return ged::Options::EditCosts::LETTER;
	}
	else {
		return ged::Options::EditCosts::CHEM_1;
	}
}

std::string dir(const std::string & dataset) {
	std::string root_dir("../../data/datasets/");
	if ((dataset == "AIDS") or (dataset == "Mutagenicity")) {
		return (root_dir + dataset + "/data/");
	}
	else if (dataset == "Letter") {
		return (root_dir + dataset + "/HIGH/");
	}
	else {
		throw ged::Error("Invalid dataset specified. Usage: ./median_tests <AIDS|Mutagenicity|Letter>");
	}
	return "";
}

std::string train_collection(const std::string & dataset) {
	std::string collection_file("../collections/");
	collection_file += dataset + "-train.xml";
	return collection_file;
}

std::string test_collection(const std::string & dataset) {
	std::string collection_file("../collections/");
	collection_file += dataset + "-test.xml";
	return collection_file;
}


int main(int argc, char* argv[]) {

	if (argc <= 1) {
		throw ged::Error("No dataset specified. Usage: ./median_chem <AIDS|Mutagenicity|Letter>");
	}
	std::string dataset(argv[1]);

	// Generate the result file.
	std::string result_filename("../output/");
	result_filename += dataset + "_INDEXING_RESULTS.csv";
	std::ofstream result_file(result_filename.c_str());
	result_file << "quantile,scan_time,query_time\n";
	result_file.close();


	// Set up the environment.
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(edit_costs(dataset));
	std::vector<ged::GEDGraph::GraphID> train_graph_ids(env.load_gxl_graphs(dir(dataset), train_collection(dataset),
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));
	std::vector<ged::GEDGraph::GraphID> test_graph_ids(env.load_gxl_graphs(dir(dataset), test_collection(dataset),
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));

	std::vector<ged::GEDGraph::GraphID> focal_graph_ids;
	for (std::size_t counter{0}; counter < 2 * train_graph_ids.size(); counter++) {
		focal_graph_ids.emplace_back(env.add_graph(dataset + "_BST_" + std::to_string(counter) + ".gxl", "no_class"));
	}
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);

	// Estimate the GEDs to define the thresholds.
	std::vector<double> estimated_ged;
	for (ged::GEDGraph::GraphID g_id : test_graph_ids) {
		for (ged::GEDGraph::GraphID h_id : test_graph_ids) {
			env.set_method(ged::Options::GEDMethod::BRANCH_FAST);
			env.run_method(g_id, h_id);
			estimated_ged.emplace_back((env.get_lower_bound(g_id, h_id) + env.get_upper_bound(g_id, h_id)) / 2);
		}
	}
	std::sort(estimated_ged.begin(), estimated_ged.end());
	std::vector<double> quantiles{0.01, 0.025, 0.05, 0.075, 0.1};
	std::vector<double> thresholds;
	for (double quantile : quantiles) {
		double real_index{quantile * static_cast<double>(estimated_ged.size() - 1)};
		std::size_t index{static_cast<std::size_t>(real_index)};
		double frac{real_index - static_cast<double>(index)};
		double threshold{estimated_ged.at(index)};
		if (index + 1 < estimated_ged.size()) {
			threshold = (1 - frac) * estimated_ged.at(index) + frac * estimated_ged.at(index + 1);
		}
		thresholds.emplace_back(threshold);
	}


	// Set up and initialize the GraphBST.
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> median_estimator(&env, constant_node_costs(dataset));
	median_estimator.set_options("--time-limit 600 --stdout 0 --init-type RANDOM --random-inits 8 --refine FALSE");
	ged::GraphBST<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> graph_bst(&env, &median_estimator);
	graph_bst.set_options("--stdout 1");
	graph_bst.init(train_graph_ids, focal_graph_ids);
	graph_bst.save("../data/indexing/" + dataset + "_BST.ini", "../data/indexing"); // Save the tree, just in case something goes wrong.
	graph_bst.set_lower_bound_method(ged::Options::GEDMethod::BLP_NO_EDGE_LABELS, "--relax TRUE --project-to-node-map FALSE");
	graph_bst.set_upper_bound_method(ged::Options::GEDMethod::IPFP);

	// Run the queries.
	for (std::size_t quantile_pos{0}; quantile_pos < 5; quantile_pos++) {
		double threshold{thresholds.at(quantile_pos)};
		std::cout << "Running tests for quantile " << quantiles.at(quantile_pos) << ".\n";
		ged::ProgressBar progress(test_graph_ids.size());
		// Run linear scan.
		std::cout << "\rRunning linear scans:" << progress << std::flush;
		double mean_scan_time{0};
		for (ged::GEDGraph::GraphID query_id : test_graph_ids) {
			auto start = std::chrono::high_resolution_clock::now();
			std::vector<ged::GEDGraph::GraphID> filtered_graphs;
			std::vector<ged::GEDGraph::GraphID> verified_graphs;
			std::vector<ged::GEDGraph::GraphID> undecided_graphs;
			for (ged::GEDGraph::GraphID data_id : train_graph_ids) {
				env.set_method(ged::Options::GEDMethod::BLP_NO_EDGE_LABELS, "--relax TRUE --project-to-node-map FALSE");
				env.run_method(query_id, data_id);
				if (env.get_lower_bound(query_id, data_id) > threshold) {
					filtered_graphs.emplace_back(data_id);
					continue;
				}
				env.set_method(ged::Options::GEDMethod::IPFP);
				env.run_method(query_id, data_id);
				if (env.get_upper_bound(query_id, data_id) <= threshold) {
					verified_graphs.emplace_back(data_id);
				}
				else {
					undecided_graphs.emplace_back(data_id);
				}
				std::cout << "\rRunning linear scans:" << progress << std::flush;
			}
			ged::Seconds runtime(std::chrono::high_resolution_clock::now() - start);
			mean_scan_time += runtime.count();
			progress.increment();
			std::cout << "\rRunning linear scans:" << progress << std::flush;
		}
		std::cout << "\n";
		mean_scan_time /= static_cast<double>(test_graph_ids.size());
		// Run query against GraphBST.
		double mean_query_time{0};
		for (ged::GEDGraph::GraphID query_id : test_graph_ids) {
			graph_bst.process_range_query(query_id, threshold);
			mean_query_time += graph_bst.get_query_time();
		}
		mean_query_time /= static_cast<double>(test_graph_ids.size());
		result_file.open(result_filename.c_str(),std::ios_base::app);
		result_file << quantiles.at(quantile_pos) << "," << mean_scan_time << "," << mean_query_time << "\n";
		result_file.close();
	}

}
