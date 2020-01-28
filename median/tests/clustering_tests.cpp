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
 * @file clustering_tests.cpp
 * @brief Tests the GraphClusteringHeuristic.
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

std::string collection(const std::string & dataset, const std::string & id) {
	std::string collection_file("../collections/");
	collection_file += dataset;
	if (dataset == "Mutagenicity") {
		collection_file += "-Correct";
	}
	return collection_file + "-90-" + id + ".xml";
}

std::size_t num_clusters(const std::string & dataset) {
	if (dataset == "Letter") {
		return 15;
	}
	return 2;
}

int main(int argc, char* argv[]) {

	if (argc <= 1) {
		throw ged::Error("No dataset specified. Usage: ./median_chem <AIDS|Mutagenicity|Letter>");
	}
	std::string dataset(argv[1]);

	// Varied sub-collections.
	std::vector<std::string> ids{"0", "1", "2", "3", "4"};

	// Generate the result file.
	std::string result_filename("../output/");
	result_filename += dataset + "_CLUSTERING_RESULTS.csv";
	std::ofstream result_file(result_filename.c_str());
	result_file << "id,method,ari,gini,sil,sod,time\n";
	result_file.close();

	ged::ProgressBar progress(ids.size() * 2);
	std::cout << "\rRunning clustering tests:" << progress << std::flush;

	for (const auto & id : ids) {

		// Set up the environment.
		ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
		env.set_edit_costs(edit_costs(dataset));
		std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(dir(dataset), collection(dataset, id),
				ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));
		std::vector<ged::GEDGraph::GraphID> focal_graph_ids;
		for (std::size_t counter{0}; counter < num_clusters(dataset); counter++) {
			focal_graph_ids.emplace_back(env.add_graph("median_" + std::to_string(counter) + ".gxl"));
		}
		env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);

		// Compute the ground truth clustering.
		std::map<std::string, std::vector<ged::GEDGraph::GraphID>> ground_truth_map;
		for (ged::GEDGraph::GraphID graph_id : graph_ids) {
			auto it = ground_truth_map.find(env.get_graph_class(graph_id));
			if (it != ground_truth_map.end()) {
				it->second.emplace_back(graph_id);
			}
			else {
				ground_truth_map.emplace(env.get_graph_class(graph_id), std::vector<ged::GEDGraph::GraphID>{graph_id});
			}
		}
		std::vector<std::vector<ged::GEDGraph::GraphID>> ground_truth_clustering;
		for (const auto & key_val : ground_truth_map) {
			ground_truth_clustering.emplace_back(key_val.second);
		}

		// Set up the clustering heuristic.
		ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> median_estimator(&env, constant_node_costs(dataset));
		median_estimator.set_options("--time-limit 600 --stdout 0 --init-type RANDOM --random-inits 8 --refine FALSE");
		ged::GraphClusteringHeuristic<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> clustering_heuristic(&env, &median_estimator);
		clustering_heuristic.set_main_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads 6");

		// Run the tests.
		std::vector<std::string> methods{"RANDOM", "MEDIAN"};
		for (const auto & method : methods) {
			if (method == "RANDOM") {
				clustering_heuristic.set_options("--refine FALSE --stdout 0 --max-itrs 0 --random-inits 1");
			}
			else {
				clustering_heuristic.set_options("--refine FALSE --stdout 0 --max-itrs 10 --random-inits 1");
			}
			clustering_heuristic.run(graph_ids, focal_graph_ids);
			result_file.open(result_filename.c_str(),std::ios_base::app);
			result_file << id << "," << method << ",";
			result_file << clustering_heuristic.get_adjusted_rand_index(ground_truth_clustering) << ",";
			result_file << clustering_heuristic.get_gini_coefficient() << ",";
			result_file << clustering_heuristic.get_silhouette_score() << ",";
			result_file << clustering_heuristic.get_sum_of_distances() << ",";
			result_file << clustering_heuristic.get_runtime() << "\n";
			result_file.close();
			progress.increment();
			std::cout << "\rRunning clustering tests:" << progress << std::flush;
		}
	}
	std::cout << "\n";

}
