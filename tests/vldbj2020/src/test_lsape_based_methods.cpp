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
 * @file tests/vldbj2019/src/test_lsape_based_methods.cpp
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
	std::size_t num_solutions_;
	std::string centralities_;
	std::string set_distances_;
	std::string ml_method_;


	std::string options_(const std::string & dataset) const {
		std::string options("");
		options += "--threads 6 --max-num-solutions " + std::to_string(num_solutions_) + " --centrality-method " + centralities_;
		if (set_distances_ != "") {
			options += " --led-method " + set_distances_;
			options += " --load ../ini/" + dataset + "_ring_" + set_distances_ + ".ini";
		}
		if (ml_method_ != "") {
			std::string ged_method_name{(ged_method_ == ged::Options::GEDMethod::RING_ML) ? "_ring_ml_" : "_bipartite_ml_BIPARTITE_"};
			options += " --ml-method " + ml_method_;
			if (ml_method_ == "DNN") {
				options += " --load ../ini/" + dataset + ged_method_name + "dnn.ini";
			}
			else {
				options += " --load ../ini/" + dataset + ged_method_name + "one_class_svm.ini";
			}
		}
		if (ged_method_ == ged::Options::GEDMethod::WALKS) {
			options += " --load ../ini/" + dataset + "_walks.ini";
		}
		if (ged_method_ == ged::Options::GEDMethod::SUBGRAPH) {
			options += " --load ../ini/" + dataset + "_subgraph.ini";
		}
		return options;
	}

public:
	Method(ged::Options::GEDMethod ged_method, std::size_t num_solutions, std::string centralities, const std::string & set_distances = "", const std::string ml_method = "") :
		ged_method_{ged_method},
		num_solutions_{num_solutions},
		centralities_{centralities},
		set_distances_{set_distances},
		ml_method_{ml_method} {
			if (num_solutions_ <= 0) {
				throw ged::Error("Invalid number of solutions.");
			}
			if (centralities_ != "NONE" and centralities_ != "DEGREE" and centralities_ != "EIGENVALUE" and centralities_ != "PAGERANK") {
				throw ged::Error("Invalid node centralities.");
			}
			if (ged_method_ == ged::Options::GEDMethod::RING) {
				if (set_distances != "LSAPE_OPTIMAL" and set_distances != "LSAPE_GREEDY" and set_distances != "GAMMA") {
					throw ged::Error("Invalid set distances.");
				}
			}
			else if (set_distances_ != "") {
				throw ged::Error("Invalid set distances.");
			}
			if (ged_method_ == ged::Options::GEDMethod::BIPARTITE_ML or ged_method_ == ged::Options::GEDMethod::RING_ML) {
				if (ml_method != "DNN" and ml_method != "ONE_CLASS_SVM") {
					throw ged::Error("Invalid machine learning method.");
				}
			}
			else if (ml_method_ != "") {
				throw ged::Error("Invalid machine learning method.");
			}
		}

		std::string name() const {
			std::stringstream name;
			if (ged_method_ == ged::Options::GEDMethod::RING) {
				if (set_distances_ == "GAMMA") {
					name << "RINGMS";
				}
				else {
					name << "RINGOPT";
				}
			}
			else if (ged_method_ == ged::Options::GEDMethod::RING_ML) {
				if (ml_method_ == "DNN") {
					name << "RINGMLDNN";
				}
				else {
					name << "RINGMLSVM";
				}
			}
			else if (ged_method_ == ged::Options::GEDMethod::BIPARTITE_ML) {
				if (ml_method_ == "DNN") {
					name << "PREDICTDNN";
				}
				else {
					name << "PREDICTSVM";
				}
			}
			else if (ged_method_ == ged::Options::GEDMethod::BIPARTITE) {
				name << "BP";
			}
			else if (ged_method_ == ged::Options::GEDMethod::BRANCH_FAST) {
				name << "BRANCHFAST";
			}
			else if (ged_method_ == ged::Options::GEDMethod::BRANCH) {
				name << "BRANCH";
			}
			else if (ged_method_ == ged::Options::GEDMethod::BRANCH_UNIFORM) {
				name << "BRANCHUNI";
			}
			else if (ged_method_ == ged::Options::GEDMethod::STAR) {
				name << "STAR";
			}
			else if (ged_method_ == ged::Options::GEDMethod::SUBGRAPH) {
				name << "SUBGRAPH";
			}
			else if (ged_method_ == ged::Options::GEDMethod::NODE) {
				name << "NODE";
			}
			else if (ged_method_ == ged::Options::GEDMethod::WALKS) {
				name << "WALKS";
			}
			name << ",$(" << num_solutions_;
			if (centralities_ == "NONE") {
				name << ",0.0)$";
			}
			else {
				name << ",0.7)$";
			}
			return name.str();
		}

		void run_on_dataset(const std::string & dataset, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env, double & avg_lb, double & avg_ub, double & avg_runtime,
				double & classification_coefficient_lb, double & classification_coefficient_ub) const {
			env.set_method(ged_method_, options_(dataset));
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


void test_on_dataset(const std::string & dataset, bool only_lb) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	util::setup_environment(dataset, false, env);

	// Collect all tested methods.
	std::vector<ged::Options::GEDMethod> ged_methods;
	if (only_lb) {
		ged_methods = {ged::Options::GEDMethod::BRANCH_UNIFORM, ged::Options::GEDMethod::NODE, ged::Options::GEDMethod::BRANCH_FAST, ged::Options::GEDMethod::BRANCH};
	}
	else {
		ged_methods = {ged::Options::GEDMethod::BIPARTITE, ged::Options::GEDMethod::STAR, ged::Options::GEDMethod::BRANCH_UNIFORM, ged::Options::GEDMethod::NODE, ged::Options::GEDMethod::WALKS, ged::Options::GEDMethod::BRANCH_FAST, ged::Options::GEDMethod::BRANCH, ged::Options::GEDMethod::SUBGRAPH, ged::Options::GEDMethod::BIPARTITE_ML, ged::Options::GEDMethod::RING, ged::Options::GEDMethod::RING_ML};
	}
	std::vector<std::string> set_distances{"GAMMA", "LSAPE_OPTIMAL"};
	std::vector<std::string> centralities{"NONE", "PAGERANK"};
	std::vector<std::string> ml_methods{"DNN", "ONE_CLASS_SVM"};
	std::vector<std::size_t> nums_solutions{1, 4, 7, 10};
	std::vector<Method> methods;
	for (auto ged_method : ged_methods) {
		for (const auto & centrality_method : centralities) {
			for (auto num_solutions : nums_solutions) {
				if (ged_method == ged::Options::GEDMethod::RING) {
					for (const auto & set_distance : set_distances) {
						methods.emplace_back(ged_method, num_solutions, centrality_method, set_distance);
					}
				}
				else if (ged_method == ged::Options::GEDMethod::RING_ML or ged_method == ged::Options::GEDMethod::BIPARTITE_ML) {
					for (const auto & ml_method : ml_methods) {
						methods.emplace_back(ged_method, num_solutions, centrality_method, "", ml_method);
					}
				}
				else {
					methods.emplace_back(ged_method, num_solutions, centrality_method);
				}
			}
		}
	}

	// Run the tests.
	std::string result_filename("../results/");
	if (only_lb) {
		result_filename += dataset + "__lsape_based_methods_LB.csv";
	}
	else {
		result_filename += dataset + "__lsape_based_methods.csv";
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
	std::vector<std::string> datasets;
	bool only_lb{false};
	int i{1};
	if (argc > 1) {
		if (std::string(argv[i]) == "--lb") {
			only_lb = true;
			i++;
		}
	}
	for (; i < argc; i++) {
		datasets.push_back(std::string(argv[i]));
		util::check_dataset(datasets.back());
	}
	if (datasets.empty()) {
		util::setup_datasets(datasets);
	}
	for (auto dataset : datasets) {
		try {
			test_on_dataset(dataset, only_lb);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << ".\n";
		}
	}
	return 0;
}



