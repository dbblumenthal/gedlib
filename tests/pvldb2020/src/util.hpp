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


#ifndef TESTS_PVLDB2020_SRC_UTIL_HPP_
#define TESTS_PVLDB2020_SRC_UTIL_HPP_

#define GXL_GEDLIB_SHARED
#include "../../../src/env/ged_env.hpp"

namespace util {


bool is_chemical_dataset(const std::string & dataset) {
	return ((dataset == "AIDS") or (dataset == "Mutagenicity") or (dataset == "acyclic") or (dataset == "alkane") or (dataset == "mao") or (dataset == "pah") );
}

bool is_letter_dataset(const std::string & dataset) {
	return ((dataset == "Letter_HIGH") or (dataset == "Letter_LOW") or (dataset == "Letter_MED"));
}

void check_dataset(const std::string & dataset) {
	if (not (is_chemical_dataset(dataset) or is_letter_dataset(dataset) or (dataset == "CMU-GED") or (dataset == "Fingerprint") or (dataset == "GREC") or (dataset == "Protein"))) {
		throw ged::Error(std::string("Dataset \"") + dataset + "\" does not exists.");
	}
}

std::string graph_dir(const std::string & dataset) {
	std::string root_dir("../../../data/datasets/");
	if ((dataset == "AIDS") or (dataset == "Fingerprint") or (dataset == "GREC") or (dataset == "Protein") or (dataset == "Mutagenicity")) {
		return (root_dir + dataset + "/data/");
	}
	else if ((dataset == "Letter_HIGH")) {
		return (root_dir + "Letter/HIGH/");
	}
	else if ((dataset == "Letter_LOW")) {
		return (root_dir + "Letter/LOW/");
	}
	else if ((dataset == "Letter_MED")) {
		return (root_dir + "Letter/MED/");
	}
	else if (dataset == "CMU-GED") {
		return (root_dir + dataset + "/CMU/");
	}
	else if ((dataset == "acyclic") or (dataset == "alkane") or (dataset == "mao") or (dataset == "pah")) {
		return (root_dir + dataset + "/");
	}
	else {
		throw ged::Error(std::string("Dataset \"") + dataset + "\" does not exists.");
	}
	return "";
}

std::string collection(const std::string & dataset) {
    std::string root_dir("../../../data/collections/");
    check_dataset(dataset);
    if (is_letter_dataset(dataset)) {
        return root_dir + "Letter.xml";
    }
    if (dataset == "Mutagenicity") {
        return root_dir + "Mutagenicity-Correct.xml";
    }
    return root_dir + dataset + ".xml";
}

std::string config_prefix(const std::string & dataset) {
	check_dataset(dataset);
	return std::string("../output/" + dataset + "_");
}

ged::Options::EditCosts edit_costs(const std::string & dataset) {
	if (is_chemical_dataset(dataset)) {
		return ged::Options::EditCosts::CHEM_2;
	}
	else if (is_letter_dataset(dataset)) {
		return ged::Options::EditCosts::LETTER;
	}
	else if (dataset == "CMU-GED") {
		return ged::Options::EditCosts::CMU;
	}
	else if (dataset == "Fingerprint") {
		return ged::Options::EditCosts::FINGERPRINT;
	}
	else if (dataset == "GREC") {
		return ged::Options::EditCosts::GREC_2;
	}
	else if (dataset == "Protein") {
		return ged::Options::EditCosts::PROTEIN;
	}
	else {
		throw ged::Error(std::string("Dataset \"") + dataset + "\" does not exists.");
	}
	return ged::Options::EditCosts::CONSTANT;
}

ged::Options::GXLNodeEdgeType node_type(const std::string & dataset) {
	check_dataset(dataset);
	if ((dataset == "Fingerprint") or (dataset == "CMU-GED")) {
		return ged::Options::GXLNodeEdgeType::UNLABELED;
	}
	return ged::Options::GXLNodeEdgeType::LABELED;
}

ged::Options::GXLNodeEdgeType edge_type(const std::string & dataset) {
	check_dataset(dataset);
	if (is_letter_dataset(dataset)) {
		return ged::Options::GXLNodeEdgeType::UNLABELED;
	}
	return ged::Options::GXLNodeEdgeType::LABELED;
}

std::unordered_set<std::string> irrelevant_node_attributes(const std::string & dataset) {
	check_dataset(dataset);
	std::unordered_set<std::string> irrelevant_attributes;
	if ((dataset == "AIDS")) {
		irrelevant_attributes.insert({"x", "y", "symbol", "charge"});
	}
	else if (dataset == "Protein") {
		irrelevant_attributes.insert("aaLength");
	}
	return irrelevant_attributes;
}

std::unordered_set<std::string> irrelevant_edge_attributes(const std::string & dataset) {
	check_dataset(dataset);
	std::unordered_set<std::string> irrelevant_attributes;
	if ((dataset == "GREC")) {
		irrelevant_attributes.insert({"angle0", "angle1"});
	}
	else if (dataset == "Protein") {
		irrelevant_attributes.insert({"distance0", "distance1"});
	}
	else if (dataset == "Fingerprint") {
		irrelevant_attributes.insert("angle");
	}
	return irrelevant_attributes;
}

ged::Options::InitType init_type(const std::string & dataset) {
	if (is_chemical_dataset(dataset) or (dataset == "Protein")) {
		return ged::Options::InitType::EAGER_WITH_SHUFFLED_COPIES;
	}
	return ged::Options::InitType::LAZY_WITH_SHUFFLED_COPIES;
}

std::vector<ged::GEDGraph::GraphID> setup_environment(const std::string & dataset, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env) {
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir(dataset), collection(dataset), node_type(dataset), edge_type(dataset), irrelevant_node_attributes(dataset), irrelevant_edge_attributes(dataset)));
	env.set_edit_costs(edit_costs(dataset));
	env.init(init_type(dataset));
	return graph_ids;
}

void setup_datasets(std::vector<std::string> & datasets) {
	datasets = {"Letter_HIGH", "pah", "AIDS", "Protein", "GREC", "Fingerprint"};
}

}

#endif /* TESTS_PVLDB2020_SRC_UTIL_HPP_ */
