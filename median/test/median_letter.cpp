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
 * @file median_letter.cpp
 * @brief Computes median graphs for Letter graphs.
 */
#define GXL_GEDLIB_SHARED

#include "../src/median_graph_estimator.hpp"

void save_letter_graph_as_tikz_file(const ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & median, const std::string & tikz_file_name) {
	std::ofstream tikz_file(tikz_file_name.c_str());
	tikz_file << "\\documentclass[tikz]{standalone}\n";
	tikz_file << "\\begin{document}\n";
	tikz_file << "\\begin{tikzpicture}\n";
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		tikz_file << "\\node[anchor=base,fill,circle,inner sep=1pt] (" << i << ") at (" << median.node_labels.at(i).at("x") << "," << median.node_labels.at(i).at("y") << ") {};\n";
	}
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		for (std::size_t j{i + 1}; j < median.num_nodes; j++) {
			if (median.adj_matrix[i][j] == 1) {
				tikz_file << "\\draw (" << i << ") -- (" << j << ");\n";
			}
		}
	}
	tikz_file << "\\end{tikzpicture}\n";
	tikz_file << "\\end{document}\n";
	tikz_file.close();
}

int main(int argc, char* argv[]) {

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(ged::Options::EditCosts::LETTER);
	std::string letter_class("A");
	if (argc > 1) {
		letter_class = std::string(argv[1]);
	}
	std::string seed("0");
	if (argc > 2) {
		seed = std::string(argv[2]);
	}
	std::string collection_file("../collections/Letter_" + letter_class + ".xml");
	std::string graph_dir("../../data/datasets/Letter/HIGH/");
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir, collection_file,
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::UNLABELED));
	ged::GEDGraph::GraphID median_id{env.add_graph("median_Letter_HIGH_" + letter_class, letter_class)};
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
	std::string ipfp_options("--threads 6 --initial-solutions 5 --initialization-method RANDOM");
	env.set_method(ged::Options::GEDMethod::IPFP, ipfp_options);
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> median_estimator(&env, false);
	median_estimator.set_options("--init-type RANDOM --randomness PSEUDO --seed " + seed);
	median_estimator.run(graph_ids, median_id);
	std::string gxl_file_name("../output/gen_median_Letter_HIGH_" + letter_class + ".gxl");
	env.save_as_gxl_graph(median_id, gxl_file_name);
	std::string tikz_file_name("../output/gen_median_Letter_HIGH_" + letter_class + ".tex");
	save_letter_graph_as_tikz_file(env.get_graph(median_id), tikz_file_name);
}



