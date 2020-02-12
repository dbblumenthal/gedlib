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
 * @file  tune_ibd_costs.cpp
 * @brief Tunes IBDCosts.
 */

#define GXL_GEDLIB_SHARED
#include "../../../src/env/ged_env.hpp"
#include "ibd_costs.hpp"
#include "CLI11.hpp"

int main(int argc, char* argv[]) {

	// Parse command line argument.
	CLI::App app;
	ged::Options::GEDMethod method{ged::Options::GEDMethod::BRANCH};
	std::map<std::string, ged::Options::GEDMethod> map{{"BRANCH", ged::Options::GEDMethod::BRANCH}, {"BRANCH_FAST", ged::Options::GEDMethod::BRANCH_FAST}};
	app.add_option("-m,--method", method, "Employed GED method.")->transform(CLI::CheckedTransformer(map, CLI::ignore_case));
	CLI11_PARSE(app, argc, argv);

	// Load the graphs and initialize the environment.
	std::cout << "Initializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	std::string graph_collection("../data/IBD.xml");
	std::string graph_dir("../data/");
	std::unordered_set<std::string> irrelevant_node_attributes;
	irrelevant_node_attributes.insert("Bacteria");
	IBDCosts<ged::GXLLabel, ged::GXLLabel> ibd_costs("../data/phylogenetic_OTU_distances.csv");
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir, graph_collection, ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes));
	env.set_edit_costs(&ibd_costs);
	env.init(ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES);
	env.set_method(method, "--threads 10");
	env.init_method();

	// Construct the folds for cross-validation.
	std::cout << "Constructing folds for cross-validation ...\n";
	std::vector<ged::GEDGraph::GraphID> control_ids;
	std::vector<ged::GEDGraph::GraphID> case_ids;
	for (ged::GEDGraph::GraphID graph_id : graph_ids) {
		if (env.get_graph_class(graph_id) == "0") {
			control_ids.emplace_back(graph_id);
		}
		else {
			case_ids.emplace_back(graph_id);
		}
	}
	std::mt19937 urng;
	std::shuffle(control_ids.begin(), control_ids.end(), urng);
	std::shuffle(case_ids.begin(), case_ids.end(), urng);
	std::vector<std::size_t> test_fold(env.num_graphs(), ged::undefined());
	for (std::size_t pos{0}; pos < control_ids.size(); pos++) {
		test_fold.at(control_ids.at(pos)) = pos % 5;
	}
	for (std::size_t pos{0}; pos < case_ids.size(); pos++) {
		test_fold.at(case_ids.at(pos)) = pos % 5;
	}
	for (std::size_t fold_id{0}; fold_id < 5; fold_id++) {
		std::cout << "Fold " << fold_id << ": ";
		std::size_t num_controls{0};
		std::size_t num_cases{0};
		for (ged::GEDGraph::GraphID graph_id : graph_ids) {
			if (test_fold.at(graph_id) != fold_id) {
				if (env.get_graph_class(graph_id) == "0") {
					num_controls++;
				}
				else {
					num_cases++;
				}
			}
		}
		std::cout << num_controls << " controls, " << num_cases << " cases.\n";
	}

	// Run cross-validation.
	std::cout << "Running cross-validation ...\n";
	std::vector<ged::GEDGraph::GraphID> selected_graph_ids;
	std::string filename;
	if (method == ged::Options::GEDMethod::BRANCH) {
		filename = "../output/BRANCH_TUNING_results.csv";
	}
	else {
		filename = "../output/BRANCH_FAST_TUNING_results.csv";
	}
	std::ofstream tuning_results(filename.c_str());
	tuning_results << "alpha,dunn-index\n";
	tuning_results.close();

	for (int exponent{-1}; exponent >= -5; exponent--) {
		double dunn_index{0.0};
		double alpha{1 - std::pow(2, static_cast<double>(exponent))};
		ibd_costs.set_alpha(alpha);
		std::cout << "Running cross-validation for alpha = " << alpha << ".\n";
		for (std::size_t fold_id{0}; fold_id < 5; fold_id++) {
			selected_graph_ids.clear();
			for (ged::GEDGraph::GraphID graph_id : graph_ids) {
				if (test_fold.at(graph_id) != fold_id) {
					selected_graph_ids.emplace_back(graph_id);
				}
			}
			double control_distance{0};
			double case_distance{0};
			double inter_class_distance{0};
			std::size_t num_control_pairs{0};
			std::size_t num_case_pairs{0};
			std::size_t num_inter_class_pairs{0};
			ged::ProgressBar progress((selected_graph_ids.size() * (selected_graph_ids.size() - 1)) / 2);
			std::cout << "\rFold " << fold_id << ": " << progress << std::flush;
			for (std::size_t pos_1{0}; pos_1 < selected_graph_ids.size() - 1; pos_1++) {
				for (std::size_t pos_2{pos_1 + 1}; pos_2 < selected_graph_ids.size(); pos_2++) {
					env.run_method(selected_graph_ids.at(pos_1), selected_graph_ids.at(pos_2));
					double ged{env.get_upper_bound(selected_graph_ids.at(pos_1), selected_graph_ids.at(pos_2))};
					if (env.get_graph_class(selected_graph_ids.at(pos_1)) == env.get_graph_class(selected_graph_ids.at(pos_2))) {
						if (env.get_graph_class(selected_graph_ids.at(pos_1)) == "0") {
							control_distance += ged;
							num_control_pairs++;
						}
						else {
							case_distance += ged;
							num_case_pairs++;
						}
					}
					else {
						inter_class_distance += ged;
						num_inter_class_pairs++;
					}
					progress.increment();
					std::cout << "\rFold " << fold_id << ": " << progress << std::flush;
				}
			}
			std::cout << "\n";
			control_distance /= static_cast<double>(num_control_pairs);
			case_distance /= static_cast<double>(num_case_pairs);
			inter_class_distance /= static_cast<double>(num_inter_class_pairs);
			dunn_index += inter_class_distance / std::max(control_distance, case_distance);
		}
		tuning_results.open(filename.c_str(), std::ios_base::app);
		tuning_results << alpha << "," << dunn_index / 5.0 << "\n";
		tuning_results.close();
	}

}


