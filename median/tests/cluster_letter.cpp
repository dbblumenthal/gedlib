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
 * @file cluster_letter.cpp
 * @brief
 */

#define GXL_GEDLIB_SHARED

#include "../src/graph_clustering_heuristic.hpp"


int main(int argc, char* argv[]) {

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(ged::Options::EditCosts::LETTER);
	std::string graph_dir("../../data/datasets/Letter/HIGH/");
	// std::vector<std::string> letter_classes = {"A", "E", "F", "H", "I", "K", "L", "M", "N", "T", "V", "W", "X", "Y", "Z"};
	std::vector<std::string> letter_classes = {"A"};
	std::string seed("0");
	if (argc > 1) {
		seed = std::string(argv[1]);
	}
	std::vector<std::vector<ged::GEDGraph::GraphID>> ground_truth_clustering;
	ged::ProgressBar progress(letter_classes.size());
	std::cout << "\rLoading GXL graphs: " << progress << std::flush;
	for (const auto & letter_class : letter_classes) {
		std::string collection_file("../collections/Letter_" + letter_class + ".xml");
		ground_truth_clustering.emplace_back(env.load_gxl_graphs(graph_dir, collection_file, ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::UNLABELED));
		progress.increment();
		std::cout << "\rLoading GXL graphs: " << progress << std::flush;
	}
	std::cout << "\n";
	std::vector<ged::GEDGraph::GraphID> graph_ids;
	for (const auto & cluster : ground_truth_clustering) {
		for (ged::GEDGraph::GraphID graph_id : cluster) {
			graph_ids.emplace_back(graph_id);
		}
	}
	std::vector<ged::GEDGraph::GraphID> focal_graph_ids;
	for (std::size_t class_id{0}; class_id < letter_classes.size(); class_id++) {
		focal_graph_ids.emplace_back(env.add_graph("focal_graph_" + std::to_string(class_id), "no_class"));
	}
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> median_estimator(&env, false);
	median_estimator.set_options("--stdout 0 --refine FALSE --seed " + seed);
	ged::GraphClusteringHeuristic<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> clustering_heuristic(&env, &median_estimator);
	clustering_heuristic.set_options("--focal-graphs MEDIANS --init-type CLUSTERS --random-inits 3 --seed " + seed);
	clustering_heuristic.run(graph_ids, focal_graph_ids);
	std::cout << "ARI: " << clustering_heuristic.get_adjusted_rand_index(ground_truth_clustering) << "\n";
}

