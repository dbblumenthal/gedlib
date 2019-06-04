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
 * @file tests/vldbj2019/src/test_best_methods.cpp
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
		std::string options("--threads 6");
		if (ged_method_ == ged::Options::GEDMethod::F2 or ged_method_ == ged::Options::GEDMethod::BLP_NO_EDGE_LABELS) {
			options += " --relax TRUE";
		}
		return options;
	}

public:
	Method(ged::Options::GEDMethod ged_method) :
		ged_method_{ged_method} {}

		std::string name() const {
			std::stringstream name;
			if (ged_method_ == ged::Options::GEDMethod::BRANCH_UNIFORM) {
				name << "BRANCHUNI";
			}
			else if (ged_method_ == ged::Options::GEDMethod::BRANCH_FAST) {
				name << "BRANCHFAST";
			}
			else if (ged_method_ == ged::Options::GEDMethod::NODE) {
				name << "NODE";
			}
			else if (ged_method_ == ged::Options::GEDMethod::F2) {
				name << "FTWO";
			}
			else if (ged_method_ == ged::Options::GEDMethod::BLP_NO_EDGE_LABELS) {
				name << "JUSTICEIP";
			}
			else if (ged_method_ == ged::Options::GEDMethod::REFINE) {
				name << "REFINE";
			}
			else {
				name << "IPFP";
			}
			return name.str();
		}

		void run_on_dataset(const std::string & dataset, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env, double & avg_lb, double & avg_ub, double & avg_runtime) const {
			env.set_method(ged_method_, options_());
			if (dataset != "Protein" or ged_method_ != ged::Options::GEDMethod::PARTITION) {
				env.init_method();
			}
			std::size_t num_runs{env.graph_ids().second * env.graph_ids().second};
			ged::ProgressBar progress_bar(num_runs);
			std::cout << "\r\t" << name() << ": " << progress_bar << std::flush;
			avg_runtime = 0;
			avg_ub = 0;
			avg_lb = 0;
			for (ged::GEDGraph::GraphID g_id = env.graph_ids().first; g_id != env.graph_ids().second; g_id++) {
				for (ged::GEDGraph::GraphID h_id = env.graph_ids().first; h_id != env.graph_ids().second; h_id++) {
					env.run_method(g_id, h_id);
					avg_lb += env.get_lower_bound(g_id, h_id);
					avg_ub += env.get_upper_bound(g_id, h_id);
					avg_runtime += env.get_runtime(g_id, h_id);

					progress_bar.increment();
					std::cout << "\r\t" << name() << ": " << progress_bar << std::flush;
				}
			}
			avg_lb /= static_cast<double>(num_runs);
			avg_ub /= static_cast<double>(num_runs);
			avg_runtime /= static_cast<double>(num_runs);
			std::cout << "\n";
		}
};


void test_on_dataset(const std::string & dataset) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";

	// Collect all tested methods.
	std::vector<ged::Options::GEDMethod> ged_methods{ged::Options::GEDMethod::NODE, ged::Options::GEDMethod::BRANCH_UNIFORM, ged::Options::GEDMethod::BRANCH_FAST, ged::Options::GEDMethod::IPFP, ged::Options::GEDMethod::REFINE, ged::Options::GEDMethod::F2, ged::Options::GEDMethod::BLP_NO_EDGE_LABELS};
	std::vector<Method> methods;
	for (auto ged_method : ged_methods) {
		methods.emplace_back(ged_method);
	}

	// Collect all suffixes.
	std::size_t max_max_size_div_10{0};
	if (dataset == "AIDS") {
		max_max_size_div_10 = 8;
	}
	else if (dataset == "Protein") {
		max_max_size_div_10 = 6;
	}
	else if (dataset == "Mutagenicity") {
		max_max_size_div_10 = 10;
	}

	// Write the header of the result file.
	std::string result_filename("../results/");
	result_filename += dataset + "__best_methods.csv";
	std::ofstream result_file(result_filename.c_str());
	for (const auto & method : methods) {
		result_file << method.name() + "_avg_lb," + method.name() + "_avg_ub," + method.name() + "_avg_runtime,";
	}
	result_file << "avg_num_nodes\n";
	result_file.close();
	// Run the tests.
	for (std::size_t max_size_dev_10{1}; max_size_dev_10 <= max_max_size_div_10; max_size_dev_10++) {
		ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
		util::setup_environment(dataset, max_size_dev_10, env);
		double avg_ub{0};
		double avg_lb{0};
		double avg_runtime{0};
		for (auto & method : methods) {
			method.run_on_dataset(dataset, env, avg_lb, avg_ub, avg_runtime);
			result_file.open(result_filename.c_str(),std::ios_base::app);
			result_file << avg_lb << "," << avg_ub << "," << avg_runtime << ",";
			result_file.close();
		}
		result_file.open(result_filename.c_str(),std::ios_base::app);
		result_file << env.get_avg_num_nodes() << "\n";
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
		util::setup_size_test_datasets(datasets);
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



