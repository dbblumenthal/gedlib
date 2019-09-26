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

	// Set dataset identifiers.
	std::vector<std::string> percents{"10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};
	std::vector<std::string> ids{"0", "1", "2", "3", "4"};

	// Set options of MGE.
	std::vector<std::string> init_types{" MEDOID", " MIN", " MAX", " MEAN", " RANDOM --random-inits 1", " RANDOM --random-inits 2", " RANDOM --random-inits 4", " RANDOM --random-inits 8", " RANDOM --random-inits 16", " RANDOM --random-inits 32"};

	// Set algorithms.
	std::vector<ged::Options::GEDMethod> algos{ged::Options::GEDMethod::BRANCH_FAST, ged::Options::GEDMethod::REFINE, ged::Options::GEDMethod::REFINE, ged::Options::GEDMethod::IPFP, ged::Options::GEDMethod::IPFP};
	std::vector<std::string> algo_options{"", "", " --initial-solutions 10 --ratio-runs-from-initial-solutions .5", "", " --initial-solutions 10 --ratio-runs-from-initial-solutions .5"};

	for (const auto & percent : percents) {
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
				mge.set_options("--time-limit 600 --init-type" + init_type);
				for (std::size_t algo_id{0}; algo_id < algos.size(); algo_id++) {
					mge.set_init_method(algos.at(algo_id), "--threads 6" + algo_options.at(algo_id));
					mge.set_descent_method(algos.at(algo_id), "--threads 6" + algo_options.at(algo_id));
				}
			}
		}
	}
}
