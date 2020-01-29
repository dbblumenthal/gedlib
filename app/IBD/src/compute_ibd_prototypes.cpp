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
 * @file  compute_ibd_prototypes.cpp
 * @brief Computes prototypical graphs for both patients with and without IBD.
 */

#define GXL_GEDLIB_SHARED
#include "../../../median/src/median_graph_estimator.hpp"
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
	std::cout << "Initializing the environment and the estimator ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	std::string graph_collection("../data/IBD.xml");
	std::string graph_dir("../data/");
	std::unordered_set<std::string> irrelevant_node_attributes;
	irrelevant_node_attributes.insert("Bacteria");
	IBDCosts<ged::GXLLabel, ged::GXLLabel> ibd_costs("../data/phylogenetic_OTU_distances.csv");
	std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir, graph_collection, ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_node_attributes));
	env.set_edit_costs(&ibd_costs);
	ged::GEDGraph::GraphID control_median_id{env.add_graph("CONTROL_median.gxl", "0")};
	ged::GEDGraph::GraphID case_median_id{env.add_graph("CASE_median.gxl", "1")};
	env.init(ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES);
	ged::MedianGraphEstimator<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> mge(&env, false);
	mge.set_options("--init-type RANDOM --random-inits 8 --refine FALSE");
	mge.set_descent_method(method, "--threads 10");

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
	std::cout << "number of controls: " << control_ids.size() << ", number of cases: " << case_ids.size() << std::endl;
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
	std::vector<ged::GEDGraph::GraphID> train_control_ids;
	std::vector<ged::GEDGraph::GraphID> train_case_ids;
	std::vector<ged::GEDGraph::GraphID> test_ids;
	std::ofstream validation_results("../output/MEDIAN_validation_results.csv");
	validation_results << "fold,control_prototype,case_prototype,precision,recall,selectivity,balanced_accuracy\n";
	validation_results.close();
	for (std::size_t fold_id{0}; fold_id < 5; fold_id++) {
		train_control_ids.clear();
		train_case_ids.clear();
		test_ids.clear();
		for (ged::GEDGraph::GraphID graph_id : graph_ids) {
			if (test_fold.at(graph_id) == fold_id) {
				test_ids.emplace_back(graph_id);
			}
			else if (env.get_graph_class(graph_id) == "0") {
				train_control_ids.emplace_back(graph_id);
			}
			else {
				train_case_ids.emplace_back(graph_id);
			}
		}
		std::cout << "number of test graphs: " << test_ids.size() << ", number of train controls: " << train_control_ids.size() << ", number of train cases: " << train_case_ids.size() << std::endl;
		mge.run(train_control_ids, control_median_id);
		mge.run(train_case_ids, case_median_id);
		std::size_t true_positive{0};
		std::size_t false_positive{0};
		std::size_t true_negative{0};
		std::size_t false_negative{0};
		env.set_method(method, "--threads 10");
		for (ged::GEDGraph::GraphID graph_id : test_ids) {
			env.run_method(control_median_id, graph_id);
			env.run_method(case_median_id, graph_id);
			if (env.get_graph_class(graph_id) == "0") {
				if (env.get_upper_bound(control_median_id, graph_id) < env.get_upper_bound(case_median_id, graph_id)) {
					true_negative++;
				}
				else {
					false_positive++;
				}
			}
			else {
				if (env.get_upper_bound(control_median_id, graph_id) > env.get_upper_bound(case_median_id, graph_id)) {
					true_positive++;
				}
				else {
					false_negative++;
				}
			}
		}
		double precision{static_cast<double>(true_positive) / static_cast<double>(true_positive + false_positive)};
		double recall{static_cast<double>(true_positive) / static_cast<double>(true_positive + false_negative)};
		double selectivity{static_cast<double>(true_negative) / static_cast<double>(true_negative + false_positive)};
		std::string control_gxl_file(std::to_string(fold_id) + "_" + env.get_graph_name(control_median_id));
		std::string case_gxl_file(std::to_string(fold_id) + "_" + env.get_graph_name(case_median_id));
		validation_results.open("../output/MEDIAN_validation_results.csv", std::ios_base::app);
		validation_results << fold_id << "," << control_gxl_file << "," << case_gxl_file << ",";
		validation_results << precision << "," << recall << "," << selectivity << "," << (recall + selectivity) / 2.0 << "\n";
		validation_results.close();
		env.save_as_gxl_graph(control_median_id, std::string("../output/") + control_gxl_file);
		env.save_as_gxl_graph(case_median_id, std::string("../output/") + case_gxl_file);
	}

	// Compute representatives for all graphs contained in the classes.
	std::cout << "Computing overall prototypes ...\n";
	mge.run(control_ids, control_median_id);
	mge.run(case_ids, case_median_id);
	env.save_as_gxl_graph(control_median_id, std::string("../output/") + env.get_graph_name(control_median_id));
	env.save_as_gxl_graph(case_median_id, std::string("../output/") + env.get_graph_name(case_median_id));

}
