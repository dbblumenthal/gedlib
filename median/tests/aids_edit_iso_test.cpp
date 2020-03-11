/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2020 by David B. Blumenthal                              *
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
 * @file  aids_edit_iso_test.cpp
 * @brief 
 */

#define GXL_GEDLIB_SHARED

#include "../../src/env/ged_env.hpp"

int main(int argc, char* argv[]) {
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(ged::Options::EditCosts::CONSTANT);
	std::string dataset("../../data/datasets/AIDS-EDIT/");
	std::string collection("../../data/collections/AIDS-EDIT.xml");
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(dataset, collection));
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
	env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM, "--threads 6");
	env.init_method();
	ged::IMatrix is_isomorphic(graph_ids.size(), graph_ids.size(), 1);
	ged::ProgressBar progress((graph_ids.size() * (graph_ids.size() - 1)) / 2);
	std::cout << "\rChecking isomorphism " << progress << std::flush;
	for (std::size_t row{0}; row < graph_ids.size() - 1; row++) {
		for (std::size_t col{row + 1}; col < graph_ids.size(); col++) {
			env.run_method(graph_ids.at(row), graph_ids.at(col));
			if (env.get_lower_bound(graph_ids.at(row), graph_ids.at(col)) > 0) {
				is_isomorphic(row, col) = 0;
				is_isomorphic(col, row) = 0;
			}
			progress.increment();
			std::cout << "\rChecking isomorphism " << progress << std::flush;
		}
	}
	std::cout << "\n";
	std::string result_filename("../output/AIDS-EDIT-ISO.csv");
	std::ofstream result_file(result_filename.c_str());
	for (auto graph_id : graph_ids) {
		result_file << "," << env.get_graph_name(graph_id);
	}
	result_file << "\n";
	for (std::size_t row{0}; row < graph_ids.size(); row++) {
		result_file << env.get_graph_name(graph_ids.at(row));
		for (std::size_t col{0}; col < graph_ids.size(); col++) {
			result_file << "," << is_isomorphic(row, col);
		}
		result_file << "\n";
	}
	result_file.close();
}

