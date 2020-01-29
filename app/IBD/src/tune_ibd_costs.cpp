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

	// Run cross-validation.
	std::cout << "Running cross-validation ...\n";
	std::vector<ged::GEDGraph::GraphID> data_ids;
	std::vector<ged::GEDGraph::GraphID> test_ids;
	std::string filename;
	if (method == ged::Options::GEDMethod::BRANCH) {
		filename = "../output/BRANCH_TUNING_results.csv";
	}
	else {
		filename = "../output/BRANCH_FAST_TUNING_results.csv";
	}
	std::ofstream tuning_results(filename.c_str());
	tuning_results << "alpha,precision,recall,selectivity,balanced_accuracy\n";
	tuning_results.close();

	for (int exponent{-1}; exponent >= -5; exponent--) {
		double alpha{1 - std::pow(2, static_cast<double>(exponent))};
		ibd_costs.set_alpha(alpha);
		double precision{0};
		double recall{0};
		double selectivity{0};
		double balanced_accuracy{0};
		std::cout << "Running cross-validation for alpha = " << alpha << ".\n";
		for (std::size_t fold_id{0}; fold_id < 5; fold_id++) {
			data_ids.clear();
			test_ids.clear();
			for (ged::GEDGraph::GraphID graph_id : graph_ids) {
				if (test_fold.at(graph_id) == fold_id) {
					test_ids.emplace_back(graph_id);
				}
				else {
					data_ids.emplace_back(graph_id);
				}
			}
			std::size_t true_positive{0};
			std::size_t false_positive{0};
			std::size_t true_negative{0};
			std::size_t false_negative{0};
			ged::ProgressBar progress(test_ids.size() * data_ids.size());
			std::cout << "\rFold " << fold_id << ": " << progress << std::flush;
			for (ged::GEDGraph::GraphID graph_id : test_ids) {
				ged::GEDGraph::GraphID closest_data_graph_id{ged::undefined()};
				double distance_to_closests_data_graph{std::numeric_limits<double>::infinity()};
				for (ged::GEDGraph::GraphID data_graph_id : data_ids) {
					env.run_method(data_graph_id, graph_id);
					progress.increment();
					std::cout << "\rFold " << fold_id << ": " << progress << std::flush;
					if (env.get_upper_bound(data_graph_id, graph_id) < distance_to_closests_data_graph) {
						closest_data_graph_id = data_graph_id;
						distance_to_closests_data_graph = env.get_upper_bound(data_graph_id, graph_id);
					}
				}
				if (env.get_graph_class(graph_id) == "0") {
					if (env.get_graph_class(closest_data_graph_id) == "0") {
						true_negative++;
					}
					else {
						false_positive++;
					}
				}
				else {
					if (env.get_graph_class(closest_data_graph_id) == "1") {
						true_positive++;
					}
					else {
						false_negative++;
					}
				}
			}
			std::cout << "\n";
			precision += static_cast<double>(true_positive) / static_cast<double>(true_positive + false_positive);
			recall += static_cast<double>(true_positive) / static_cast<double>(true_positive + false_negative);
			selectivity += static_cast<double>(true_negative) / static_cast<double>(true_negative + false_positive);
			balanced_accuracy += (static_cast<double>(true_positive) / static_cast<double>(true_positive + false_negative) + static_cast<double>(true_negative) / static_cast<double>(true_negative + false_positive)) / 2.0;
		}
		tuning_results.open(filename.c_str(), std::ios_base::app);
		tuning_results << alpha << "," << precision / 5.0 << "," << recall / 5.0 << "," << selectivity / 5.0 << "," << balanced_accuracy / 5.0 << "\n";
		tuning_results.close();
	}

}


