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
 * @file tests/ijprai2020/src/util.hpp
 * @brief Provides utility functions for tests of PR submission.
 * @details This file was used for the experiments in the following submission:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun:
 *   &ldquo;Upper Bounding GED via Transformations to LSAPE Based on Rings and Machine Learning&rdquo;,
 *   Submitted to TKDE.
 */

#ifndef SRC_TESTS_PR2018_UTIL_HPP_
#define SRC_TESTS_PR2018_UTIL_HPP_

#define GXL_GEDLIB_SHARED
#include "../../../src/env/ged_env.hpp"

namespace util {

bool is_synth_mol_dataset(const std::string & dataset) {
	return ((dataset == "S-MOL_NL01") or (dataset == "S-MOL_NL04") or (dataset == "S-MOL_NL07") or (dataset == "S-MOL_NL10"));
}

bool is_synth_acyclic_dataset(const std::string & dataset) {
    return ((dataset == "S-acyclic_NL03") or (dataset == "S-acyclic_NL05") or (dataset == "S-acyclic_NL07") or (dataset == "S-acyclic_NL09"));
}

bool is_synth_mao_dataset(const std::string & dataset) {
	return ((dataset == "S-mao_NL03") or (dataset == "S-mao_NL05") or (dataset == "S-mao_NL07") or (dataset == "S-mao_NL09"));
}

bool is_chemical_dataset(const std::string & dataset) {
	return (is_synth_mao_dataset(dataset) or is_synth_mol_dataset(dataset) or is_synth_acyclic_dataset(dataset) or (dataset == "AIDS") or (dataset == "Mutagenicity") or (dataset == "acyclic") or (dataset == "alkane") or (dataset == "mao") or (dataset == "pah") );
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
	else if (dataset == "S-MOL_NL01") {
		return (root_dir + "S-MOL/NL01/");
	}
	else if (dataset == "S-MOL_NL04") {
		return (root_dir + "S-MOL/NL04/");
	}
	else if (dataset == "S-MOL_NL07") {
		return (root_dir + "S-MOL/NL07/");
	}
	else if (dataset == "S-MOL_NL10") {
		return (root_dir + "S-MOL/NL10/");
	}
	else if (dataset == "S-mao_NL03") {
		return (root_dir + "S-mao/NL03/");
	}
	else if (dataset == "S-mao_NL05") {
		return (root_dir + "S-mao/NL05/");
	}
	else if (dataset == "S-mao_NL07") {
		return (root_dir + "S-mao/NL07/");
	}
	else if (dataset == "S-mao_NL09") {
		return (root_dir + "S-mao/NL09/");
	}
	else if (dataset == "S-acyclic_NL03") {
		return (root_dir + "S-acyclic/NL03/");
	}
	else if (dataset == "S-acyclic_NL05") {
		return (root_dir + "S-acyclic/NL05/");
	}
	else if (dataset == "S-acyclic_NL07") {
		return (root_dir + "S-acyclic/NL07/");
	}
	else if (dataset == "S-acyclic_NL09") {
		return (root_dir + "S-acyclic/NL09/");
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

std::string train_collection(const std::string & dataset) {
	std::string root_dir("../collections/");
	check_dataset(dataset);
	if (is_letter_dataset(dataset)) {
		return (root_dir + "Letter_50.xml");
	}
	if (is_synth_mol_dataset(dataset)) {
		return (root_dir + "S-MOL_50.xml");
	}
	if (is_synth_acyclic_dataset(dataset)) {
		return (root_dir + "acyclic_50.xml");
	}
	if (is_synth_mao_dataset(dataset)) {
		return (root_dir + "mao_50.xml");
	}
	return root_dir + dataset + "_50.xml";
}

std::string test_collection(const std::string & dataset) {
	std::string root_dir("../collections/");
	check_dataset(dataset);
	if (is_letter_dataset(dataset)) {
		return (root_dir + "Letter_100.xml");
	}
	if (is_synth_mol_dataset(dataset)) {
		return (root_dir + "S-MOL_100.xml");
	}
	if (is_synth_acyclic_dataset(dataset)) {
		return (root_dir + "acyclic_100.xml");
	}
	if (is_synth_mao_dataset(dataset)) {
		return (root_dir + "mao_100.xml");
	}
	return root_dir + dataset + "_100.xml";
}

std::string config_prefix(const std::string & dataset) {
	check_dataset(dataset);
	return std::string("../output/" + dataset + "_");
}

std::string init_options(const std::string & dataset, const std::string & config_suffix, const std::string & data_suffix = "", bool save_train = false, bool load_train = false, std::size_t threads = 10) {
	check_dataset(dataset);
	std::string options("--threads ");
	options += std::to_string(threads) + " --save ../output/";
	options += dataset + "_" + config_suffix + ".ini";
	if (save_train) {
		if (load_train) {
			throw ged::Error("Training data cannot be both saved and loaded.");
		}
		options += " --save-train ../output/" + dataset + "_" + data_suffix + ".data";
	}
	if (load_train) {
		options += " --load-train ../output/" + dataset + "_" + data_suffix + ".data";
	}
	return options;
}

std::string ground_truth_option(const std::string & dataset) {
	check_dataset(dataset);
	return std::string(" --ground-truth-method IPFP");
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
		return ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES;
	}
	return ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES;
}

std::vector<ged::GEDGraph::GraphID> setup_environment(const std::string & dataset, bool train, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env) {
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir(dataset), (train ? train_collection(dataset) : test_collection(dataset)), node_type(dataset), edge_type(dataset), irrelevant_node_attributes(dataset), irrelevant_edge_attributes(dataset)));
	env.set_edit_costs(edit_costs(dataset));
	env.init(init_type(dataset));
	return graph_ids;
}

void setup_datasets(std::vector<std::string> & datasets) {
	datasets = {"Letter_HIGH", "pah", "AIDS", "Protein", "GREC", "Fingerprint"};
}

}

#endif /* SRC_TESTS_PR2018_UTIL_HPP_ */
