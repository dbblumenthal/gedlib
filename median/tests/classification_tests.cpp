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

#include "../src/graph_clustering_heuristic.hpp"

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
	result_filename += dataset + "_CLASSIFICATION_RESULTS.csv";
	std::ofstream result_file(result_filename.c_str());
	result_file << "method,classif,time\n";
	result_file.close();


	// Set up the environment.
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(edit_costs(dataset));
	std::vector<ged::GEDGraph::GraphID> train_graph_ids(env.load_gxl_graphs(dir(dataset), train_collection(dataset),
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));
	std::vector<ged::GEDGraph::GraphID> test_graph_ids(env.load_gxl_graphs(dir(dataset), test_collection(dataset),
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
	env.set_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads 6");

	double mean_runtime{0};
	double mean_classif_coeff{0};
	for (auto h_id : test_graph_ids) {
		ged::GEDGraph::GraphID closest_graph_id{ged::undefined()};
		double dist_to_closest_graph{std::numeric_limits<double>::infinity()};
		auto start = std::chrono::high_resolution_clock::now();
		for (auto g_id : train_graph_ids) {
			env.run_method(g_id, h_id);
			if (env.get_upper_bound(g_id, h_id) < dist_to_closest_graph) {
				dist_to_closest_graph = env.get_upper_bound(g_id, h_id);
				closest_graph_id = g_id;
			}
		}
		ged::Seconds runtime(std::chrono::high_resolution_clock::now() - start);
		mean_runtime += runtime.count();
		if (env.get_graph_class(h_id) == env.get_graph_class(closest_graph_id)) {
			mean_classif_coeff += 1;
		}
	}
	mean_runtime /= static_cast<double>(test_graph_ids.size());
	mean_classif_coeff /= static_cast<double>(test_graph_ids.size());
	result_file.open(result_filename.c_str(),std::ios_base::app);
	result_file << "SCAN," << mean_classif_coeff << "," << mean_runtime << "\n";
	result_file.close();

	// Compute the training classes.
	std::map<std::string, std::vector<ged::GEDGraph::GraphID>> training_classes;
	for (ged::GEDGraph::GraphID graph_id : train_graph_ids) {
		auto it = training_classes.find(env.get_graph_class(graph_id));
		if (it != training_classes.end()) {
			it->second.emplace_back(graph_id);
		}
		else {
			training_classes.emplace(env.get_graph_class(graph_id), std::vector<ged::GEDGraph::GraphID>{graph_id});
		}
	}

	// Set up the clustering heuristic.
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> median_estimator(&env, constant_node_costs(dataset));
	median_estimator.set_options("--time-limit 600 --stdout 0 --init-type RANDOM --random-inits 8 --refine FALSE");
	ged::GraphClusteringHeuristic<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> clustering_heuristic(&env, &median_estimator);
	clustering_heuristic.set_main_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads 6");

	// Run the tests.
	ged::ProgressBar progress(training_classes.size());
	std::cout << "\rRunning classification tests:" << progress << std::flush;
	std::vector<ged::GEDGraph::GraphID> all_focal_graph_ids;
	for (const auto & key_val : training_classes) {
		std::vector<ged::GEDGraph::GraphID> focal_graph_ids;
		for (std::size_t counter{0}; counter < 10; counter++) {
			focal_graph_ids.emplace_back(env.add_graph(dataset + "_" + key_val.first + "_" + std::to_string(counter) + ".gxl", key_val.first));
			all_focal_graph_ids.emplace_back(focal_graph_ids.back());
		}
		env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
		clustering_heuristic.set_options("--refine FALSE --stdout 0 --max-itrs 10 --random-inits 1");
		clustering_heuristic.run(key_val.second, focal_graph_ids);
		progress.increment();
		std::cout << "\rRunning classification tests:" << progress << std::flush;
	}

	mean_runtime = 0;
	mean_classif_coeff = 0;
	for (auto h_id : test_graph_ids) {
		ged::GEDGraph::GraphID closest_graph_id{ged::undefined()};
		double dist_to_closest_graph{std::numeric_limits<double>::infinity()};
		auto start = std::chrono::high_resolution_clock::now();
		for (auto g_id : all_focal_graph_ids) {
			env.run_method(g_id, h_id);
			if (env.get_upper_bound(g_id, h_id) < dist_to_closest_graph) {
				dist_to_closest_graph = env.get_upper_bound(g_id, h_id);
				closest_graph_id = g_id;
			}
		}
		ged::Seconds runtime(std::chrono::high_resolution_clock::now() - start);
		mean_runtime += runtime.count();
		if (env.get_graph_class(h_id) == env.get_graph_class(closest_graph_id)) {
			mean_classif_coeff += 1;
		}
	}
	mean_runtime /= static_cast<double>(test_graph_ids.size());
	mean_classif_coeff /= static_cast<double>(test_graph_ids.size());
	result_file.open(result_filename.c_str(),std::ios_base::app);
	result_file << "MEDIAN" << "," << mean_classif_coeff << "," << mean_runtime << "\n";
	result_file.close();
	std::cout << "\n";

}
