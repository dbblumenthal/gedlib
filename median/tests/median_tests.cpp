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
 * @file median_tests.cpp
 * @brief Tests the MedianGraphEstimator.
 */

#define GXL_GEDLIB_SHARED

#include "../src/median_graph_estimator.hpp"

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

std::string collection(const std::string & dataset, const std::string & percent, const std::string & id) {
	std::string collection_file("../collections/");
	collection_file += dataset;
	if (dataset == "Mutagenicity") {
		collection_file += "-Correct";
	}

	return collection_file + "-" + percent + "-" + id + ".xml";
}

int main(int argc, char* argv[]) {

	if (argc <= 1) {
		throw ged::Error("No dataset specified. Usage: ./median_chem <AIDS|Mutagenicity|Letter>");
	}
	std::string dataset(argv[1]);

	// Varied sub-collections.
	std::vector<std::string> percents{"10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};
	std::vector<std::string> ids{"4", "0", "1", "2", "3"};

	// Varied estimator parameters.
	std::vector<std::string> init_types{"RANDOM", "MAX", "MIN", "MEAN", "MEDOID"};
	std::vector<std::string> nums_inits{"16", "1", "2", "4", "8", "32"};

	// Varied algorithm parameters.
	std::vector<ged::Options::GEDMethod> algos{ged::Options::GEDMethod::IPFP, ged::Options::GEDMethod::BRANCH_FAST, ged::Options::GEDMethod::REFINE};
	std::vector<std::string> algo_options_suffixes{" --initial-solutions 10 --ratio-runs-from-initial-solutions .5", "", " --initial-solutions 10 --ratio-runs-from-initial-solutions .5"};

	// Generate the result file.
	std::string result_filename("../output/");
	result_filename += dataset + "_RESULTS.csv";
	std::ofstream result_file(result_filename.c_str());
	result_file << "percent,id,init_type,num_inits,algo,time,time_init,time_converged,sod,sod_init,sod_converged,itrs,state\n";
	result_file.close();

	// Iterate through all varied sub-collections.
	for (const auto & percent : percents) {

		// Initialize progress bar.
		ged::ProgressBar progress(ids.size() * (init_types.size() + nums_inits.size() - 1) * algos.size());
		std::cout << "\rRunning tests for " << percent << "% collections:" << progress << std::flush;

		for (const auto & id : ids) {

			// Set up the environment.
			ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
			env.set_edit_costs(edit_costs(dataset));
			std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(dir(dataset), collection(dataset, percent, id),
					ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));
			ged::GEDGraph::GraphID median_id{env.add_graph("median")};
			env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);

			// Set up the estimator.
			ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> mge(&env, constant_node_costs(dataset));
			mge.set_refine_method(ged::Options::GEDMethod::IPFP, "--threads 6 --initial-solutions 10 --ratio-runs-from-initial-solutions .5");

			for (const auto & init_type : init_types) {
				for (const auto & num_inits : nums_inits) {
					std::string mge_options("--time-limit 600 --stdout 0 --init-type " + init_type);
					if (init_type != "RANDOM" and num_inits != "1") {
						continue;
					}
					else {
						std::random_device rng;
						mge_options += " --random-inits " + num_inits + " --seed " + std::to_string(rng());
					}
					for (std::size_t algo_id{0}; algo_id < algos.size(); algo_id++) {
						// Select the GED algorithm.
						ged::Options::GEDMethod algo{algos.at(algo_id)};
						std::string algo_options("--threads 6" + algo_options_suffixes.at(algo_id));
						mge.set_options(mge_options);
						mge.set_init_method(algo, algo_options);
						mge.set_descent_method(algo, algo_options);

						// Run the estimator.
						mge.run(graph_ids, median_id);

						// Write the results.
						result_file.open(result_filename.c_str(),std::ios_base::app);
						result_file << percent << "," << id << "," << init_type << "," << num_inits << "," << algo;
						result_file << "," << mge.get_runtime() << "," << mge.get_runtime(ged::Options::AlgorithmState::INITIALIZED) << "," << mge.get_runtime(ged::Options::AlgorithmState::CONVERGED);
						result_file << "," << mge.get_sum_of_distances() << "," << mge.get_sum_of_distances(ged::Options::AlgorithmState::INITIALIZED) << "," << mge.get_sum_of_distances(ged::Options::AlgorithmState::CONVERGED);
						std::vector<std::size_t> nums_itrs(mge.get_num_itrs());
						result_file << "," << nums_itrs.at(0);
						for (std::size_t pos{1}; pos < nums_itrs.size(); pos++) {
							result_file << ";" << nums_itrs.at(pos);
						}
						result_file << "," << mge.get_state() << "\n";
						result_file.close();

						// Increment the progress bar and print current progress.
						progress.increment();
						std::cout << "\rRunning tests for " << percent << "% collections:" << progress << std::flush;
					}
				}
			}
		}
		std::cout << "\n";
	}
}
