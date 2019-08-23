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
 * @file tests/vldbj2019/src/test_ls_based_methods.cpp
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
	std::size_t num_initial_solutions_;
	double ratio_initial_solutions_;
	double randpost_penalty_;
	std::size_t num_randpost_loops_;
	std::size_t max_swap_size_;
	std::size_t num_orderings_;


	std::string options_() const {
		std::string options("");
		options += "--threads 6 --initial-solutions " + std::to_string(num_initial_solutions_) + " --ratio-runs-from-initial-solutions " + std::to_string(ratio_initial_solutions_) + " --num-randpost-loops " + std::to_string(num_randpost_loops_) + " --randpost-penalty " + std::to_string(randpost_penalty_);
		if (ged_method_ == ged::Options::GEDMethod::REFINE) {
			options += " --max-swap-size " + std::to_string(max_swap_size_);
		}
		else if (ged_method_ == ged::Options::GEDMethod::BP_BEAM) {
			options += " --num-orderings " + std::to_string(num_orderings_);
		}
		return options;
	}

public:
	Method(ged::Options::GEDMethod ged_method, std::size_t num_initial_solutions, double ratio_initial_solutions, std::size_t num_randpost_loops, double randpost_penalty = 0.0, std::size_t max_swap_size = 2, std::size_t num_orderings = 1) :
		ged_method_{ged_method},
		num_initial_solutions_{num_initial_solutions},
		ratio_initial_solutions_{ratio_initial_solutions},
		randpost_penalty_{randpost_penalty},
		num_randpost_loops_{num_randpost_loops},
		max_swap_size_{max_swap_size},
		num_orderings_{num_orderings} {}

		std::string name() const {
			std::stringstream name;
			if (ged_method_ == ged::Options::GEDMethod::BP_BEAM) {
				if (num_orderings_ == 1) {
					name << "BPBEAM";
				}
				else {
					name << "IBPBEAM";
				}
			}
			else if (ged_method_ == ged::Options::GEDMethod::REFINE) {
				if (max_swap_size_ == 2) {
					name << "REFINE";
				}
				else {
					name << "KREFINE";
				}
			}
			else {
				name << "IPFP";
			}
			name << ",$(" << num_initial_solutions_ << "," << ratio_initial_solutions_ << "," << num_randpost_loops_ << "," << randpost_penalty_ << ")$";
			return name.str();
		}

		void run_on_dataset(const std::string & dataset, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env, double & avg_lb, double & avg_ub, double & avg_runtime,
				double & classification_coefficient_lb, double & classification_coefficient_ub) const {
			env.set_method(ged_method_, options_());
			env.init_method();
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

struct RandpostSetup {
	std::size_t num_initial_solutions;
	double ratio_initial_solutions;
	std::size_t num_randpost_loops;
	RandpostSetup (std::size_t num_initial_solutions, double ratio_initial_solutions, std::size_t num_randpost_loops) :
		num_initial_solutions{num_initial_solutions},
		ratio_initial_solutions{ratio_initial_solutions},
		num_randpost_loops{num_randpost_loops} {}
};


void test_on_dataset(const std::string & dataset, bool only_refine, bool only_ipfp, bool only_bp_beam) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	util::setup_environment(dataset, false, env);

	// Collect all tested methods.
	std::vector<ged::Options::GEDMethod> ged_methods;
	if (only_refine) {
		ged_methods = {ged::Options::GEDMethod::REFINE};
	}
	else if (only_ipfp) {
		ged_methods = {ged::Options::GEDMethod::IPFP};
	}
	else if (only_bp_beam) {
		ged_methods = {ged::Options::GEDMethod::BP_BEAM};
	}
	else {
		ged_methods = {ged::Options::GEDMethod::REFINE, ged::Options::GEDMethod::IPFP, ged::Options::GEDMethod::BP_BEAM};
	}
	std::vector<std::size_t> nums_orderings{1, 20};
	std::vector<std::size_t> max_swap_sizes{2, 3};
	std::vector<RandpostSetup> randpost_setups;
	randpost_setups.emplace_back(1, 1, 0);
	randpost_setups.emplace_back(10, 1, 0);
	randpost_setups.emplace_back(20, 1, 0);
	randpost_setups.emplace_back(30, 1, 0);
	randpost_setups.emplace_back(40, 1, 0);
	randpost_setups.emplace_back(40, 0.5, 1);
	randpost_setups.emplace_back(40, 0.25, 3);
	randpost_setups.emplace_back(40, 0.125, 7);
	std::vector<double> randpost_penalties{0.0, 1.0};

	std::vector<Method> methods;
	for (auto ged_method : ged_methods) {
		for (auto randpost_setup : randpost_setups) {
			if (randpost_setup.num_randpost_loops > 0) {
				for (auto randpost_penalty : randpost_penalties) {
					if (ged_method == ged::Options::GEDMethod::REFINE) {
						for (auto max_swap_size : max_swap_sizes) {
							methods.emplace_back(ged_method, randpost_setup.num_initial_solutions, randpost_setup.ratio_initial_solutions, randpost_setup.num_randpost_loops, randpost_penalty, max_swap_size);
						}
					}
					else if (ged_method == ged::Options::GEDMethod::BP_BEAM) {
						for (auto num_orderings : nums_orderings) {
							methods.emplace_back(ged_method, randpost_setup.num_initial_solutions, randpost_setup.ratio_initial_solutions, randpost_setup.num_randpost_loops, randpost_penalty, 2, num_orderings);
						}
					}
					else {
						methods.emplace_back(ged_method, randpost_setup.num_initial_solutions, randpost_setup.ratio_initial_solutions, randpost_setup.num_randpost_loops, randpost_penalty);
					}
				}
			}
			else {
				if (ged_method == ged::Options::GEDMethod::REFINE) {
					for (auto max_swap_size : max_swap_sizes) {
						methods.emplace_back(ged_method, randpost_setup.num_initial_solutions, randpost_setup.ratio_initial_solutions, randpost_setup.num_randpost_loops, 0.0, max_swap_size);
					}
				}
				else if (ged_method == ged::Options::GEDMethod::BP_BEAM) {
					for (auto num_orderings : nums_orderings) {
						methods.emplace_back(ged_method, randpost_setup.num_initial_solutions, randpost_setup.ratio_initial_solutions, randpost_setup.num_randpost_loops, 0.0, 2, num_orderings);
					}
				}
				else {
					methods.emplace_back(ged_method, randpost_setup.num_initial_solutions, randpost_setup.ratio_initial_solutions, randpost_setup.num_randpost_loops);
				}
			}
		}
	}


	std::string result_filename("../results/");
	if (only_refine) {
		result_filename += dataset + "__ls_based_methods_REFINE.csv";
	}
	else if (only_ipfp) {
		result_filename += dataset + "__ls_based_methods_IPFP.csv";
	}
	else if (only_bp_beam) {
		result_filename += dataset + "__ls_based_methods_BP_BEAM.csv";
	}
	else {
		result_filename += dataset + "__ls_based_methods.csv";
	}
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
	int i{1};
	bool only_refine{false};
	bool only_ipfp{false};
	bool only_bp_beam{false};
	if (argc > 1) {
		if (std::string(argv[i]) == "--refine") {
			only_refine = true;
			i++;
		}
		else if (std::string(argv[i]) == "--ipfp") {
			only_ipfp = true;
			i++;
		}
		else if (std::string(argv[i]) == "--bp-beam") {
			only_bp_beam = true;
			i++;
		}
	}
	std::vector<std::string> datasets;
	for (; i < argc; i++) {
		datasets.push_back(std::string(argv[i]));
		util::check_dataset(datasets.back());
	}
	if (datasets.empty()) {
		util::setup_datasets(datasets);
	}
	for (auto dataset : datasets) {
		try {
			test_on_dataset(dataset, only_refine, only_ipfp, only_bp_beam);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << ".\n";
		}
	}
	return 0;
}



