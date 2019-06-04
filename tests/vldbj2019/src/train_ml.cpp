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
 * @file train_ml.cpp
 * @brief Trains the methods ged::RingML and ged::BipartiteML for different datasets.
 * @details The binary built from this file was used for the experiments in the following paper:
 * - D. B. Blumenthal, S. Bougleux, N. Boria, J. Gamper, L. Brun:
 *   &ldquo;Comparing heuristics for graph edit distance computation&rdquo;,
 *   Accepted for publication in VLDB J.
 */

#include "util.hpp"

void train_on_dataset(const std::string & dataset) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	util::setup_environment(dataset, true, env);

	// Train RingML and BipartiteML with DNN.
	env.set_method(ged::Options::GEDMethod::RING_ML, util::init_options(dataset, "ring_ml_dnn", "ring_ml", true, false, 16) + util::ground_truth_option(dataset) + " --ml-method DNN --save-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
	env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, util::init_options(dataset, "bipartite_ml_BIPARTITE_dnn", "bipartite_ml_BIPARTITE", true, false, 16) + " --ml-method DNN --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();

	// Train RingML and BipartiteML with OneClassSVM.
	env.set_method(ged::Options::GEDMethod::RING_ML, util::init_options(dataset, "ring_ml_one_class_svm", "ring_ml_one_class", true, false, 16) + util::ground_truth_option(dataset) + " --ml-method ONE_CLASS_SVM --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
	env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, util::init_options(dataset, "bipartite_ml_BIPARTITE_one_class_svm", "bipartite_ml_BIPARTITE_one_class", true, false, 16) + " --ml-method ONE_CLASS_SVM --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();

}

int main(int argc, char* argv[]) {
	std::vector<std::string> datasets;
	for (int i{1}; i < argc; i++) {
		datasets.push_back(std::string(argv[i]));
		util::check_dataset(datasets.back());
	}
	if (datasets.empty()) {
		util::setup_datasets(datasets);
	}
	for (auto dataset : datasets) {
		try {
			train_on_dataset(dataset);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << ".\n";
		}
	}
	return 0;
}
