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
 * @file tests/vldbj2019/src/test_misc_methods.cpp
 * @brief Runs tests for VLDB J. submission.
 * @details The binary built from this file was used for the experiments in the following paper:
 * - D. B. Blumenthal, N. Boria, J. Gamper, S. Bougleux L. Brun:
 *   &ldquo;Comparing heuristics for graph edit distance computation&rdquo;,
 *   VLDB J. 2019
 */

#include "util.hpp"

class Method {
private:
	// method and options
	ged::Options::GEDMethod ged_method_;


	std::string options_() const {
		std::string options("");
		if (ged_method_ != ged::Options::GEDMethod::BRANCH_COMPACT and ged_method_ != ged::Options::GEDMethod::PARTITION) {
			options += "--threads 6";
		}
		if (ged_method_ == ged::Options::GEDMethod::HYBRID) {
			options += " --time-limit 1";
		}
		return options;
	}

public:
	Method(ged::Options::GEDMethod ged_method) :
		ged_method_{ged_method} {}

		std::string name() const {
			std::stringstream name;
			if (ged_method_ == ged::Options::GEDMethod::BRANCH_COMPACT) {
				name << "BRANCHCOMPACT";
			}
			else if (ged_method_ == ged::Options::GEDMethod::PARTITION) {
				name << "PARTITION";
			}
			else if (ged_method_ == ged::Options::GEDMethod::SIMULATED_ANNEALING) {
				name << "SA";
			}
			else if (ged_method_ == ged::Options::GEDMethod::HYBRID) {
				name << "HYBRID";
			}
			else if (ged_method_ == ged::Options::GEDMethod::BRANCH_TIGHT) {
				name << "BRANCHTIGHT";
			}
			else if (ged_method_ == ged::Options::GEDMethod::HED) {
				name << "HED";
			}
			return name.str();
		}

		void run_on_dataset(const std::string & dataset, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env, double & avg_lb, double & avg_ub, double & avg_runtime,
				double & classification_coefficient_lb, double & classification_coefficient_ub) const {
			env.set_method(ged_method_, options_());
			if (dataset != "Protein" or ged_method_ != ged::Options::GEDMethod::PARTITION) {
				env.init_method();
			}
			std::size_t num_runs{env.graph_ids().second * env.graph_ids().second};
			ged::ProgressBar progress_bar(num_runs);
			std::cout << "\r\t" << name() << ": " << progress_bar << std::flush;
			std::size_t num_intra_class_runs{0};
			std::size_t num_inter_class_runs{0};
			double largest_lb{0.0};
			double avg_intra_class_lb{0.0};
			double avg_inter_class_lb{0.0};
			double largest_ub{0.0};
			double avg_intra_class_ub{0.0};
			double avg_inter_class_ub{0.0};
			avg_runtime = 0;
			avg_ub = 0;
			avg_lb = 0;
			for (ged::GEDGraph::GraphID g_id = env.graph_ids().first; g_id != env.graph_ids().second; g_id++) {
				for (ged::GEDGraph::GraphID h_id = env.graph_ids().first; h_id != env.graph_ids().second; h_id++) {
					env.run_method(g_id, h_id);
					avg_lb += env.get_lower_bound(g_id, h_id);
					avg_ub += env.get_upper_bound(g_id, h_id);
					avg_runtime += env.get_runtime(g_id, h_id);
					if (env.get_lower_bound(g_id, h_id) > largest_lb) {
						largest_lb = env.get_lower_bound(g_id, h_id);
					}
					if (env.get_upper_bound(g_id, h_id) > largest_ub) {
						largest_ub = env.get_upper_bound(g_id, h_id);
					}
					if (env.get_graph_class(g_id) == env.get_graph_class(h_id)) {
						avg_intra_class_lb += env.get_lower_bound(g_id, h_id);
						avg_intra_class_ub += env.get_upper_bound(g_id, h_id);
						num_intra_class_runs++;
					}
					else {
						avg_inter_class_lb += env.get_lower_bound(g_id, h_id);
						avg_inter_class_ub += env.get_upper_bound(g_id, h_id);
						num_inter_class_runs++;
					}
					progress_bar.increment();
					std::cout << "\r\t" << name() << ": " << progress_bar << std::flush;
				}
			}
			avg_lb /= static_cast<double>(num_runs);
			avg_ub /= static_cast<double>(num_runs);
			avg_runtime /= static_cast<double>(num_runs);
			avg_intra_class_lb /= static_cast<double>(num_intra_class_runs);
			avg_intra_class_ub /= static_cast<double>(num_intra_class_runs);
			avg_inter_class_lb /= static_cast<double>(num_inter_class_runs);
			avg_inter_class_ub /= static_cast<double>(num_inter_class_runs);
			if (largest_lb > 0) {
				classification_coefficient_lb = (avg_inter_class_lb - avg_intra_class_lb) / largest_lb;
			}
			else {
				classification_coefficient_lb = 0;
			}
			if (largest_ub > 0) {
				classification_coefficient_ub = (avg_inter_class_ub - avg_intra_class_ub) / largest_ub;
			}
			else {
				classification_coefficient_ub = 0;
			}
			std::cout << "\n";
		}
};


void test_on_dataset(const std::string & dataset) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	util::setup_environment(dataset, false, env);

	// Collect all tested methods.
	std::vector<ged::Options::GEDMethod> ged_methods{ged::Options::GEDMethod::HED, ged::Options::GEDMethod::BRANCH_COMPACT, ged::Options::GEDMethod::PARTITION, ged::Options::GEDMethod::SIMULATED_ANNEALING, ged::Options::GEDMethod::HYBRID, ged::Options::GEDMethod::BRANCH_TIGHT};
	std::vector<Method> methods;
	for (auto ged_method : ged_methods) {
		methods.emplace_back(ged_method);
	}

	// Run the tests.
	std::string result_filename("../results/");
	result_filename += dataset + "__misc_methods.csv";
	std::ofstream result_file(result_filename.c_str());
	result_file << "method;avg_lb;avg_ub;avg_runtime;classification_coefficient_lb;classification_coefficient_ub\n";
	result_file.close();
	double avg_ub{0};
	double avg_lb{0};
	double avg_runtime{0};
	double classification_coefficient_lb{0};
	double classification_coefficient_ub{0};
	for (auto & method : methods) {
		method.run_on_dataset(dataset, env, avg_lb, avg_ub, avg_runtime, classification_coefficient_lb, classification_coefficient_ub);
		result_file.open(result_filename.c_str(),std::ios_base::app);
		result_file << method.name() << ";" << avg_lb << ";" << avg_ub << ";" << avg_runtime << ";" << classification_coefficient_lb << ";" << classification_coefficient_ub << "\n";
		result_file.close();
	}
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
			test_on_dataset(dataset);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << ".\n";
		}
	}
	return 0;
}



