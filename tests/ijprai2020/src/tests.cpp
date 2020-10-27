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
 * @file tests/ijprai2020/src/tests.cpp
 * @brief Carries out the tests for PR submission.
 * @details The binary built from this file was used for the experiments in the following submission:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun:
 *   &ldquo;Upper Bounding GED via Transformations to LSAPE Based on Rings and Machine Learning&rdquo;,
 *   Submitted to TKDE.
 */

#include "util.hpp"

class Method {
private:
	// method and options
	ged::Options::GEDMethod ged_method_;
	std::size_t num_threads_;
	std::size_t num_solutions_;
	std::string centralities_;
	std::string set_distances_;
	std::string ml_method_;


	std::string options_(const std::string & dataset) const {
		std::string options("");
		options += "--threads " + std::to_string(num_threads_) + " --max-num-solutions " + std::to_string(num_solutions_) + " --centrality-method " + centralities_;
		if (set_distances_ != "") {
			options += " --led-method " + set_distances_;
			options += " --load ../output/" + dataset + "_ring_" + set_distances_ + ".ini";
		}
		if (ml_method_ != "") {
			std::string ged_method_name{(ged_method_ == ged::Options::GEDMethod::RING_ML) ? "_ring_ml_" : "_bipartite_ml_BIPARTITE_"};
			if (ml_method_ == "DNN" or ml_method_ == "SVM") {
				options += " --ml-method " + ml_method_;
			}
			else if (ml_method_ == "ONE_CLASS_SVM_LIKELIHOOD") {
				options += " --ml-method ONE_CLASS_SVM --one-class-svm-likelihood TRUE";
			}
			else {
				options += " --ml-method ONE_CLASS_SVM --one-class-svm-likelihood FALSE";
			}
			if (ml_method_ == "DNN") {
				options += " --load ../output/" + dataset + ged_method_name + "dnn.ini";
			}
			else if (ml_method_ == "SVM") {
				options += " --load ../output/" + dataset + ged_method_name + "svm.ini";
			}
			else {
				options += " --load ../output/" + dataset + ged_method_name + "one_class_svm.ini";
			}
		}
		if (ged_method_ == ged::Options::GEDMethod::WALKS) {
			options += " --load ../output/" + dataset + "_walks.ini";
		}
		if (ged_method_ == ged::Options::GEDMethod::SUBGRAPH) {
			options += " --load ../output/" + dataset + "_subgraph.ini";
		}
		return options;
	}

public:
	Method(ged::Options::GEDMethod ged_method, std::size_t threads, std::size_t solutions, std::string centralities, const std::string & set_distances = "", const std::string ml_method = "") :
		ged_method_{ged_method},
		num_threads_{threads},
		num_solutions_{solutions},
		centralities_{centralities},
		set_distances_{set_distances},
		ml_method_{ml_method} {
			if (num_threads_ <= 0) {
				throw ged::Error("Invalid number of threads.");
			}
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
				if (ml_method != "DNN" and ml_method != "SVM" and ml_method != "ONE_CLASS_SVM_LIKELIHOOD" and ml_method != "ONE_CLASS_SVM_SCALE") {
					throw ged::Error("Invalid machine learning method.");
				}
			}
			else if (ml_method_ != "") {
				throw ged::Error("Invalid machine learning method.");
			}
		}

		std::string name() const {
			std::stringstream name;
			name << ged_method_ << "__C-" << centralities_;
			if (set_distances_ != "") {
				name << "__LED-" << set_distances_;
			}
			if (ml_method_ != "") {
				name << "__ML-" << ml_method_;
			}
			return name.str();
		}

		std::size_t num_threads() const {
			return num_threads_;
		}

		std::size_t num_solutions() const {
			return num_solutions_;
		}

		void run_on_dataset(const std::string & dataset, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env, double & avg_ub, double & avg_runtime, double & avg_classification_ratio) const {
			env.set_method(ged_method_, options_(dataset));
			env.init_method();
			std::size_t num_runs{(env.graph_ids().second * env.graph_ids().second) - env.graph_ids().second};
			ged::ProgressBar progress_bar(num_runs);
			std::cout << "\r\t" << name() << ", " << num_threads_ << " threads, " << num_solutions_ << " solutions: " << progress_bar << std::flush;
			ged::GEDGraph::GraphID closest_graph_id{std::numeric_limits<ged::GEDGraph::GraphID>::max()};
			double distance_to_closest_graph{std::numeric_limits<double>::infinity()};
			avg_runtime = 0;
			avg_ub = 0;
			avg_classification_ratio = 0;
			for (ged::GEDGraph::GraphID g_id = env.graph_ids().first; g_id != env.graph_ids().second; g_id++) {
				closest_graph_id = std::numeric_limits<ged::GEDGraph::GraphID>::max();
				distance_to_closest_graph = std::numeric_limits<double>::infinity();
				for (ged::GEDGraph::GraphID h_id = env.graph_ids().first; h_id != env.graph_ids().second; h_id++) {
					if (g_id == h_id) {
						continue;
					}
					env.run_method(g_id, h_id);
					avg_ub += env.get_upper_bound(g_id, h_id);
					avg_runtime += env.get_runtime(g_id, h_id);
					if (env.get_upper_bound(g_id, h_id) < distance_to_closest_graph) {
						distance_to_closest_graph = env.get_upper_bound(g_id, h_id);
						closest_graph_id = h_id;
					}
					progress_bar.increment();
					std::cout << "\r\t" << name() << ", " << num_threads() << " threads, " << num_solutions() << " solutions: " << progress_bar << std::flush;
				}
				if (env.get_graph_class(g_id) == env.get_graph_class(closest_graph_id)) {
					avg_classification_ratio += 1.0;
				}
			}
			avg_ub /= static_cast<double>(num_runs);
			avg_runtime /= static_cast<double>(num_runs);
			avg_classification_ratio /= static_cast<double>(env.graph_ids().second);
			std::cout << "\n";
		}
};


void test_on_dataset(const std::string & dataset, const std::vector<string> & ml_methods, bool quick) {

	// Initialize environment.
	std::cout << "\n=== " << dataset << " ===\n";
	std::cout << "\tInitializing the environment ...\n";
	ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
	util::setup_environment(dataset, false, env);

	// Learn the parameters.
	std::vector<ged::Options::GEDMethod> ged_methods{ged::Options::GEDMethod::BIPARTITE, ged::Options::GEDMethod::NODE, ged::Options::GEDMethod::WALKS, ged::Options::GEDMethod::BRANCH_FAST, ged::Options::GEDMethod::BRANCH, ged::Options::GEDMethod::SUBGRAPH, ged::Options::GEDMethod::BIPARTITE_ML, ged::Options::GEDMethod::RING, ged::Options::GEDMethod::RING_ML};
	std::vector<std::string> set_distances{"GAMMA", "LSAPE_GREEDY", "LSAPE_OPTIMAL"};
	std::vector<std::string> centralities{"NONE", "PAGERANK"};
	std::vector<std::size_t> threads{1, 4, 7, 10};
	std::vector<std::size_t> solutions{1, 4, 7, 10};
	if (quick) {
		centralities = {"NONE"};
		threads = {10};
		solutions = {10};
	}
	std::vector<Method> methods;
	for (auto ged_method : ged_methods) {
		for (auto num_threads : threads) {
			for (auto num_solutions : solutions) {
				for (const auto & centrality_method : centralities) {
					if (ged_method == ged::Options::GEDMethod::RING) {
						for (const auto & set_distance : set_distances) {
							methods.emplace_back(ged_method, num_threads, num_solutions, centrality_method, set_distance);
						}
					}
					else if (ged_method == ged::Options::GEDMethod::RING_ML or ged_method == ged::Options::GEDMethod::BIPARTITE_ML) {
						for (const auto & ml_method : ml_methods) {
							methods.emplace_back(ged_method, num_threads, num_solutions, centrality_method, "", ml_method);
						}
					}
					else {
						methods.emplace_back(ged_method, num_threads, num_solutions, centrality_method);
					}
				}
			}
		}
	}
	std::string result_filename("../output/");
	result_filename += dataset + "__RESULTS.csv";
	std::ofstream result_file(result_filename.c_str());
	result_file << "method,num_threads,num_solutions,avg_ub,avg_runtime,avg_classification_ratio\n";
	result_file.close();
	double avg_ub{0};
	double avg_runtime{0};
	double avg_classification_ratio{0};
	for (auto & method : methods) {
		method.run_on_dataset(dataset, env, avg_ub, avg_runtime, avg_classification_ratio);
		result_file.open(result_filename.c_str(),std::ios_base::app);
		result_file << method.name() << "," << method.num_threads() << "," << method.num_solutions() << ",";
		result_file << avg_ub << "," << avg_runtime << "," << avg_classification_ratio << "\n";
		result_file.close();
	}
}

int main(int argc, char* argv[]) {
	std::vector<std::string> datasets;
	std::vector<std::string> ml_methods{"DNN", "SVM", "ONE_CLASS_SVM_LIKELIHOOD", "ONE_CLASS_SVM_SCALE"};
	bool quick{false};
	int i{1};
	if (argc > 1) {
		std::string first_option(argv[i]);
		if (first_option == "--no-svm") {
			ml_methods = {"DNN", "ONE_CLASS_SVM_LIKELIHOOD"};
			i++;
		}
		else if (first_option == "--quick") {
			ml_methods = {"DNN", "ONE_CLASS_SVM_LIKELIHOOD"};
			quick = true;
			i++;
		}
		else {
			std::cout << "first option = \"" << first_option << "\"\n";
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
			test_on_dataset(dataset, ml_methods, quick);
		}
		catch (const std::exception & error) {
			std::cerr << error.what() << ". " << "Error on " << dataset << ".\n";
		}
	}
	return 0;
}



