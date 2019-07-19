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
 * @file median_chem.cpp
 * @brief Computes median graphs for chemical graphs.
 */
#define GXL_GEDLIB_SHARED

#include "../src/median_graph_estimator.hpp"

std::unordered_set<std::string> irrelevant_node_attributes(const std::string & dataset) {
	std::unordered_set<std::string> irrelevant_attributes;
	if ((dataset == "AIDS")) {
		irrelevant_attributes.insert({"x", "y", "symbol", "charge"});
	}
	return irrelevant_attributes;
}

std::string dir(const std::string & dataset) {
	std::string root_dir("../../data/datasets/");
	if ((dataset == "AIDS") or (dataset == "Mutagenicity")) {
		return (root_dir + dataset + "/data/");
	}
	else if ((dataset == "acyclic") or (dataset == "alkane") or (dataset == "mao") or (dataset == "pah")) {
		return (root_dir + dataset + "/");
	}
	else {
		throw ged::Error("Invalid dataset specified. Usage: ./median_chem <AIDS|Mutagenicity|acyclic|alkane|mao|pah>");
	}
	return "";
}

std::string collection(const std::string & dataset) {
	std::string root_dir("../collections/");
	return root_dir + dataset + ".xml";
}

int main(int argc, char* argv[]) {

	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	env.set_edit_costs(ged::Options::EditCosts::CHEM_1);
	if (argc <= 1) {
		throw ged::Error("No dataset specified. Usage: ./median_chem <AIDS|Mutagenicity|acyclic|alkane|mao|pah>");
	}
	std::string dataset(argv[1]);
	std::string seed("0");
	if (argc > 2) {
		seed = std::string(argv[2]);
	}
	std::string num_inits("10");
	if (argc > 3) {
		num_inits = std::string(argv[4]);
	}
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(dir(dataset), collection(dataset),
			ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes(dataset)));
	ged::GEDGraph::GraphID median_id{env.add_graph("median_" + dataset)};
	env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> median_estimator(&env, true);
	//median_estimator.set_options("--init-type MEDOID");
	median_estimator.set_options("--init-type RANDOM --refine TRUE --randomness PSEUDO --seed " + seed + " --random-inits " + num_inits);
	//median_estimator.set_descent_method(ged::Options::GEDMethod::REFINE, "--initial-solutions 5 --threads 6");
	median_estimator.set_refine_method(ged::Options::GEDMethod::IPFP, "--threads 6 --initial-solutions 40 --initialization-method RANDOM --num-randpost-loops 3");
	median_estimator.run(graph_ids, median_id);
	std::string gxl_file_name("../output/gen_median_" + dataset + ".gxl");
	env.save_as_gxl_graph(median_id, gxl_file_name);
}
