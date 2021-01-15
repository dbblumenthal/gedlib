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
 * @file  ibd_tests.cpp
 * @brief Test performance of ged::MedianGraphEstimator on IBD graphs.
 */

#define GXL_GEDLIB_SHARED

#include "../src/median_graph_estimator.hpp"
#include "../../src/edit_costs/ibd.hpp"


std::string dir() {
	return "../../data/datasets/IBD/data/";
}

std::string collection(const std::string & graph_class) {
	return "../collections/IBD_" + graph_class + ".xml";
}


int main(int argc, char* argv[]) {

	// Determine output mode.
	if (argc <= 1) {
		throw ged::Error("No output mode selected. Usage: ./ibd_tests <0|1|2>");
	}
	std::string stdout(argv[1]);
	if (stdout != "0" and stdout != "1" and stdout != "2") {
		throw ged::Error("Invalid output mode selected. Usage: ./ibd_tests <0|1|2>");
	}

	// Initialize the edit costs.
	ged::IBD<ged::GXLLabel, ged::GXLLabel> ibd_costs("../../src/edit_costs/otu_distances.csv");


	// Compute medians for all three classes.
	std::vector<std::string> graph_classes{"C", "NO", "YES"};
	for (const std::string & graph_class : graph_classes) {

		// Write header of the result file.
		std::string result_filename("../output/IBD_" + graph_class + "_NEW_RESULTS.csv");
        std::string median_filename("../output/IBD_" + graph_class + "_MEDIAN.gxl");
		std::ofstream result_file(result_filename.c_str());
		result_file << "time_bcu_I,sod_bcu_I\n";
		result_file.close();

		// Initialize the environment.
		ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
		env.set_edit_costs(&ibd_costs);
		std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(dir(), collection(graph_class)));
		ged::GEDGraph::GraphID median_id{env.add_graph("median.gxl", graph_class)};
		env.init(ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES);

		// Initialize the median graph estimator.
		ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> mge(&env, false);
		mge.set_refine_method(ged::Options::GEDMethod::IPFP, "--threads 6 --initial-solutions 10 --ratio-runs-from-initial-solutions .5");
		std::random_device rng;
		mge.set_options("--init-type RANDOM --refine TRUE --stdout " + stdout + " --random-inits 8 --randomness PSEUDO --seed " + std::to_string(rng()));


		// Run the tests and save the results.
		// mge.set_descent_method(ged::Options::GEDMethod::REFINE, "--threads 6 --initial-solutions 10 --ratio-runs-from-initial-solutions .5");
		// mge.run(graph_ids, median_id);
		// double time_bcu_II{mge.get_runtime(ged::Options::AlgorithmState::TERMINATED)};
		// double sod_bcu_II{mge.get_sum_of_distances(ged::Options::AlgorithmState::TERMINATED)};
		mge.set_descent_method(ged::Options::GEDMethod::BRANCH_FAST, "--threads 6");
		mge.run(graph_ids, median_id);
		env.save_as_gxl_graph(median_id, median_filename.c_str());
		result_file.open(result_filename.c_str(),std::ios_base::app);
		result_file << mge.get_runtime(ged::Options::AlgorithmState::TERMINATED) << "," << mge.get_sum_of_distances(ged::Options::AlgorithmState::TERMINATED) << "\n";
		// result_file << time_bcu_II << "," << sod_bcu_II << ",";
		// result_file << mge.get_runtime(ged::Options::AlgorithmState::CONVERGED) << "," << mge.get_sum_of_distances(ged::Options::AlgorithmState::CONVERGED) << "\n";
		result_file.close();
	}


}
