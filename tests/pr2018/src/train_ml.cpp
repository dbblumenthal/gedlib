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
 * @file train_ring_ml.cpp
 * @brief Trains the methods ged::RingML and ged::BipartiteML for different datasets.
 * @details The binary built from this file was used for the experiments in the following submission:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun:
 *   &ldquo;Designing heuristics for graph edit distance&rdquo;,
 *   Submitted to PR.
 */

#include "util.hpp"

void train_on_dataset(const std::string & dataset) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	util::setup_environment(dataset, true, env);

	// Train BipartiteML with matrix constructed by ged::Bipartite.
	env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, util::init_options(dataset, "bipartite_ml_BIPARTITE_dnn", "bipartite_ml_BIPARTITE", true, false, 16) + " --ml-method DNN --save-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
	env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, util::init_options(dataset, "bipartite_ml_BIPARTITE_svm", "bipartite_ml_BIPARTITE", false, true, 16) + " --ml-method SVM --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
	env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, util::init_options(dataset, "bipartite_ml_BIPARTITE_one_class_svm", "bipartite_ml_BIPARTITE_one_class", true, false, 16) + " --ml-method ONE_CLASS_SVM --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();

	// Train BipartiteML with matrix constructed by ged::Ring.
	env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, util::init_options(dataset, "bipartite_ml_RING_dnn", "bipartite_ml_RING", true, false, 16) + " --ml-method DNN --lsape-method RING --lsape-options '--load " + util::config_prefix(dataset) + "ring_LSAPE_OPTIMAL.ini' --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
	env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, util::init_options(dataset, "bipartite_ml_RING_svm", "bipartite_ml_RING", false, true, 16) + " --ml-method SVM --lsape-method RING --lsape-options '--load " + util::config_prefix(dataset) + "ring_LSAPE_OPTIMAL.ini' --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
	env.set_method(ged::Options::GEDMethod::BIPARTITE_ML, util::init_options(dataset, "bipartite_ml_RING_one_class_svm", "bipartite_ml_RING_one_class", true, false, 16) + " --ml-method ONE_CLASS_SVM --lsape-method RING --lsape-options '--load " + util::config_prefix(dataset) + "ring_LSAPE_OPTIMAL.ini' --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();

	// Train RingML.
	env.set_method(ged::Options::GEDMethod::RING_ML, util::init_options(dataset, "ring_ml_dnn", "ring_ml", true, false, 16) + util::ground_truth_option(dataset) + " --ml-method DNN --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
	env.set_method(ged::Options::GEDMethod::RING_ML, util::init_options(dataset, "ring_ml_svm", "ring_ml", false, true, 16) + util::ground_truth_option(dataset) + " --ml-method SVM --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
	env.set_method(ged::Options::GEDMethod::RING_ML, util::init_options(dataset, "ring_ml_one_class_svm", "ring_ml_one_class", true, false, 16) + util::ground_truth_option(dataset) + " --ml-method ONE_CLASS_SVM --load-ground-truth " + util::config_prefix(dataset) + "ground_truth.data");
	env.init_method();
}

int main(int argc, char* argv[]) {
	std::vector<std::string> datasets;
	util::setup_datasets(datasets);
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
