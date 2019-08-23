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
 * @file tests/sspr2018/src/learn_ring_params.cpp
 * @brief Learn the parameters of the method ged::Ring for different datasets.
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

void init_rings_on_dataset(const std::string & dataset, const std::vector<std::string> & led_methods) {
	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(std::string("../../../data/datasets/") + dataset + "/", std::string("../collections/") + dataset + "_50.xml"));
	if (dataset == "GREC") {
		env.set_edit_costs(ged::Options::EditCosts::GREC_1);
	}
	else {
		env.set_edit_costs(ged::Options::EditCosts::CHEM_1);
	}
	env.init();

	// Initialize the methods.
	for (auto led_method : led_methods) {
		std::cout << "\n=== " << led_method << " ===\n";
		env.set_method(ged::Options::GEDMethod::RING, std::string("--threads 11 --led-method ") + led_method + " --init-initial-solutions 100 --save ../output/" + dataset + "_ring_" + led_method + ".ini");
		env.init_method();
	}
}


int main(int argc, char* argv[]) {
	std::vector<std::string> led_methods{"GAMMA","LSAPE_GREEDY","LSAPE_OPTIMAL"};
	std::vector<std::string> datasets{"mao","pah","alkane","acyclic"};
	for (auto dataset : datasets) {
		try {
			init_rings_on_dataset(dataset, led_methods);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << "\n";
		}
	}
	return 0;
}
