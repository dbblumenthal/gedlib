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
 * @file bst_query_letter.cpp
 * @brief
 */

#define GXL_GEDLIB_SHARED

#include "../src/graph_bst.hpp"


int main(int argc, char* argv[]) {

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(ged::Options::EditCosts::LETTER);
	std::size_t max_cluster_size{10};
	if (argc > 1) {
		max_cluster_size = std::stoul(std::string(argv[1]));
	}
	std::string seed("0");
	if (argc > 2) {
		seed = std::string(argv[2]);
	}
	std::string query_collection_file("../collections/Letter_60.xml");
	std::string graph_dir("../../data/datasets/Letter/HIGH/");
	std::vector<ged::GEDGraph::GraphID> query_graph_ids(env.load_gxl_graphs(graph_dir, query_collection_file,
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::UNLABELED));
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);

	std::vector<double> estimated_ged;
	for (ged::GEDGraph::GraphID g_id : query_graph_ids) {
		for (ged::GEDGraph::GraphID h_id : query_graph_ids) {
			env.set_method(ged::Options::GEDMethod::BRANCH_FAST);
			env.run_method(g_id, h_id);
			estimated_ged.emplace_back((env.get_lower_bound(g_id, h_id) + env.get_upper_bound(g_id, h_id)) / 2);
		}
	}
	std::sort(estimated_ged.begin(), estimated_ged.end());



	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> mge(&env, false);
	mge.set_options("--stdout 0 --seed " + seed);
	ged::GraphBST<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> graph_bst(&env, &mge);

	std::map<double, std::vector<std::pair<std::size_t, std::size_t>>> num_evals;
	std::vector<double> thresholds;
	for (double percent : std::vector<double>({0.01, 0.025, 0.05, 0.075, 0.1})) {
		double real_index{percent * static_cast<double>(estimated_ged.size() - 1)};
		std::size_t index{static_cast<std::size_t>(real_index)};
		double frac{real_index - static_cast<double>(index)};
		double threshold{estimated_ged.at(index)};
		if (index + 1 < estimated_ged.size()) {
			threshold = (1 - frac) * estimated_ged.at(index) + frac * estimated_ged.at(index + 1);
		}
		thresholds.emplace_back(threshold);
		num_evals[threshold] = std::vector<std::pair<std::size_t, std::size_t>>();
	}
	for (std::size_t max_cluster_size : {4, 8, 16, 32, 64}) {
		std::string config("../data/Letter/Letter_150_" + std::to_string(max_cluster_size) + "_BST.ini");
		std::string focal_graph_dir("../data/Letter");
		graph_bst.load(config, graph_dir, focal_graph_dir, ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::UNLABELED);
		graph_bst.set_lower_bound_method(ged::Options::GEDMethod::BLP_NO_EDGE_LABELS);
		for (double threshold : thresholds) {
			graph_bst.process_range_query(query_graph_ids.at(0), threshold);
			num_evals[threshold].emplace_back(max_cluster_size, graph_bst.get_num_ged_evals());
		}
	}
	std::cout << "====== RESULT SUMMARY =====\n";
	for (const auto & threshold_evals : num_evals) {
		std::cout << threshold_evals.first << ":";
		for (const auto & size_evals : threshold_evals.second) {
			std::cout << " (" << size_evals.first << "," << size_evals.second << ")";
		}
		std::cout << "\n";
	}
}
