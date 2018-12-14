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
 * @file  ml_based_method.ipp
 * @brief ged::MLBasedMethod class definition.
 */

#ifndef SRC_METHODS_ML_BASED_METHOD_IPP_
#define SRC_METHODS_ML_BASED_METHOD_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
~MLBasedMethod() {
	delete ground_truth_method_;
}

template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
MLBasedMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
num_features_{undefined()},
prediction_initialized_(GEDGraph::dummy_node(), GEDGraph::dummy_node()),
assignments_(),
ml_method_{DNN},
ground_truth_method_{new IPFP<UserNodeLabel, UserEdgeLabel>(this->ged_data_)},
ground_truth_options_(""),
dnn_params_(),
dnn_(),
dnn_training_data_(),
dnn_feature_vectors_(),
dnn_types_(),
svm_params_(),
svm_(),
svm_training_data_(),
svm_feature_vectors_(),
svm_types_(),
one_class_svm_use_likelihood_{true},
one_class_svm_(),
infile_(""),
outfile_(""),
logfile_(""),
training_infile_(""),
training_outfile_(""),
ground_truth_infile_(""),
ground_truth_outfile_("") {
	this->compute_lower_bound_ = false;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
predict(const GEDGraph & g, const GEDGraph & h, const NodeMap::Assignment & assignment) {

	if ((not this->initialized_) and (not initialized_for_prediction_(g, h))) {
		lsape_pre_graph_init_(true);
		ml_init_graph_(g);
		ml_init_graph_(h);
	}
	if (not initialized_for_prediction_(g, h)) {
		ml_init_feature_variables_(g, h, this->num_threads_);
		prediction_initialized_.first = g.id();
		prediction_initialized_.second = h.id();
	}

	GEDGraph::NodeID i{assignment.first};
	GEDGraph::NodeID k{assignment.second};
	std::vector<double> feature_vector;
	if ((i != GEDGraph::dummy_node()) and (k != GEDGraph::dummy_node())) {
		ml_populate_substitution_feature_vector_(g, h, i, k, feature_vector);
	}
	else if (i != GEDGraph::dummy_node()) {
		ml_populate_deletion_feature_vector_(g, i, feature_vector);
	}
	else if (k != GEDGraph::dummy_node()) {
		ml_populate_insertion_feature_vector_(h, k, feature_vector);
	}
	else {
		return 0.0;
	}
	Assignment_ assignment_(0, 0, false, feature_vector);
	return decision_value_(assignment_);
}

// === Definitions of member functions inherited from LSAPEBasedMethod. ===

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_init_() {

	ml_init_();

	// Initialize the ground truth method.
	if (compute_or_load_ground_truth_()) {
		if (ground_truth_infile_ == "") {
			ground_truth_method_->init();
			if (ground_truth_outfile_ != "") {
				std::ofstream ofs(ground_truth_outfile_);
				ofs.close();
			}
		}
	}

	// Return if the method is initialized from a configuration file.
	if (load_config_file_()) {
		return;
	}

	// Ask the derived class about the size of the feature vectors.
	num_features_ = ml_get_num_features_();

	// Initialize the training data.
	load_or_generate_training_data_();
	if (training_outfile_ != "") {
		save_training_data_();
	}

	// Run the machine learning method.
	train_();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_pre_graph_init_(bool called_at_runtime) {
	if (load_config_file_()) {
		if (ml_method_ == DNN) {
			num_features_ = dnn_.load(infile_);
		}
		else if (ml_method_ == SVM){
			num_features_ = svm_.load(infile_);
		}
		else {
			num_features_ = one_class_svm_.load(infile_, one_class_svm_use_likelihood_);
		}
		ml_init_for_num_features_();
	}
	else if (called_at_runtime){
		throw Error("You are trying to run a machine learning based method without training it. Call init_method() before calling run_method() or provide an initialization file via the option \"--init-infile <filename>\".");
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {
	assignments_.clear();
	generate_assignments_(g, h);
	double correctly_predicted{0.0};
	double correctly_predicted_as_bad{0.0};
	double correctly_predicted_as_good{0.0};
	double ground_truth_bad{0.0};
	double ground_truth_good{0.0};
	IMatrix assignment_matrix(master_problem.num_rows(), master_problem.num_cols());
#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t pos = 0; pos < assignments_.size(); pos++) {
		Assignment_ & assignment = assignments_.at(pos);
		master_problem(assignment.row_in_master(), assignment.col_in_master()) = decision_value_(assignment);
		if (std::isnan(master_problem(assignment.row_in_master(), assignment.col_in_master()))) {
			std::string error_msg("master_problem(" + std::to_string(assignment.row_in_master()) + "," + std::to_string(assignment.col_in_master()) + ")=NaN\n");
			error_msg += "assignment=" + assignment.to_string() + "\n";
			throw Error(error_msg);
		}
		if (log_prediction_ratios_()) {
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				bool predicted_as_bad{master_problem(assignment.row_in_master(), assignment.col_in_master()) > 0.5};
				if (assignment.is_good_assignment()) {
					assignment_matrix(assignment.row_in_master(), assignment.col_in_master()) = 1;
					ground_truth_good += 1.0;
					if (not predicted_as_bad) {
						correctly_predicted += 1.0;
						correctly_predicted_as_good += 1.0;
					}
				}
				else {
					assignment_matrix(assignment.row_in_master(), assignment.col_in_master()) = 0;
					ground_truth_bad += 1.0;
					if (predicted_as_bad) {
						correctly_predicted += 1.0;
						correctly_predicted_as_bad += 1.0;
					}
				}
			}
		}
	}
	if (log_prediction_ratios_()) {
		double num_assignments{static_cast<double>(assignments_.size())};
		std::ofstream logfile;
		logfile.open(logfile_, std::ios::app);
		logfile << "=====\ng_id=" << g.id() << "\nh_id=" << h.id() << "\nassignment_matrix=\n";
		logfile << assignment_matrix << "\nmaster_problem=\n";
		logfile << master_problem << "\nbaseline_accuracy=";
		logfile << ((ground_truth_bad > ground_truth_good) ? (ground_truth_bad / num_assignments) : (ground_truth_good / num_assignments)) << "\naccuracy=";
		logfile << (correctly_predicted / num_assignments) << "\nprecision=";
		logfile << (correctly_predicted_as_bad / ground_truth_bad) << "\nrecall=";
		logfile << (correctly_predicted_as_good / ground_truth_good) << "\n";
		logfile.close();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	std::string general_options("[--ml-method <arg>] [--ground-truth-method <arg>] [--ground-truth-options <arg>] [--load <arg>] [--save <arg>] [--load-train <arg>] [--save-train <arg>] [--load-ground-truth <arg>] [--save-ground-truth <arg>] [--log <arg>]");
	std::string dnn_options("[--dnn-activation <arg>] [--dnn-hidden-layers-range <arg>] [--dnn-neurons-per-layer-range <arg>]");
	std::string svm_options("[--svm-gamma-exp-range <arg>] [--svm-c-exp-range <arg>] [--one-class-svm-likelihood <arg>]");
	if (ml_valid_options_string_() == "") {
		return (general_options + " " + dnn_options + " " + svm_options);
	}
	return (ml_valid_options_string_() + " " + general_options + " " + dnn_options + " " + svm_options);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {
	prediction_initialized_.first = GEDGraph::dummy_node();
	prediction_initialized_.second = GEDGraph::dummy_node();
	ml_method_ = DNN;
	delete ground_truth_method_;
	ground_truth_options_ = "";
	ground_truth_method_ = new IPFP<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
	ground_truth_method_->set_options(std::string("--initial-solutions 80 --ratio-runs-from-initial-solutions 0.5 --lower-bound-method BRANCH_TIGHT --threads ") + std::to_string(this->num_threads_));
	dnn_params_.activation_candidates = {FANN::activation_function_enum::RELU, FANN::activation_function_enum::SIGMOID};
	dnn_params_.min_num_hidden_layers = 1;
	dnn_params_.max_num_hidden_layers = 10;
	dnn_params_.min_num_neurons_per_layer = 1;
	dnn_params_.max_num_neurons_per_layer = 20;
	svm_params_.min_gamma_exp = -3;
	svm_params_.max_gamma_exp = 4;
	svm_params_.min_c_exp = -3;
	svm_params_.max_c_exp = 4;
	one_class_svm_use_likelihood_ = true;
	infile_ = std::string("");
	outfile_ = std::string("");
	logfile_ = std::string("");
	training_infile_ = std::string("");
	training_outfile_ = std::string("");
	ground_truth_infile_ = std::string("");
	ground_truth_outfile_ = std::string("");
	num_features_ = undefined();
	ml_set_default_options_();
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_parse_option_(const std::string & option, const std::string & arg) {
	bool is_valid_option{false};
	if (option == "ml-method") {
		if (arg == "SVM") {
			ml_method_ = SVM;
		}
		else if (arg == "ONE_CLASS_SVM") {
			ml_method_ = ONE_CLASS_SVM;
		}
		else if (arg != "DNN"){
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option ml-method. Usage: options = \"[--ml-method DNN|ONE_CLASS_SVM] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "ground-truth-method") {
		if (arg == "ANCHOR_AWARE_GED") {
			delete ground_truth_method_;
			ground_truth_method_ = new AnchorAwareGED<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
#ifdef GUROBI
		else if (arg == "F1") {
			delete ground_truth_method_;
			ground_truth_method_ = new F1<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "F2") {
			delete ground_truth_method_;
			ground_truth_method_ = new F2<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "COMPACT_MIP") {
			delete ground_truth_method_;
			ground_truth_method_ = new CompactMIP<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg != "IPFP") {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option ground-truth-method. Usage: options = \"[--ground-truth-method ANCHOR_AWARE|F1|F2|COMPACT_MIP|IPFP] [...]\"");
		}
#else
		else if (arg != "IPFP") {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option ground-truth-method. Usage: options = \"[--ground-truth-method ANCHOR_AWARE|IPFP] [...]\"");
		}
#endif
		is_valid_option = true;
	}
	else if (option == "ground-truth-options") {
		ground_truth_options_ = arg;
		std::size_t bad_option_start{ground_truth_options_.find("--threads")};
		std::size_t next_option_start;
		if (bad_option_start != std::string::npos) {
			next_option_start = ground_truth_options_.find("--", bad_option_start + 1);
			if (next_option_start != std::string::npos) {
				ground_truth_options_ = ground_truth_options_.substr(0, bad_option_start) + ground_truth_options_.substr(next_option_start);
			}
			else {
				ground_truth_options_ = ground_truth_options_.substr(0, bad_option_start);
			}
		}
		is_valid_option = true;
	}
	else if (option == "load") {
		infile_ = arg;
		is_valid_option = true;
	}
	else if (option == "save") {
		outfile_ = arg;
		is_valid_option = true;
	}
	else if (option == "load-train") {
		training_infile_ = arg;
		is_valid_option = true;
	}
	else if (option == "save-train") {
		training_outfile_ = arg;
		is_valid_option = true;
	}
	else if (option == "load-ground-truth") {
		ground_truth_infile_ = arg;
		is_valid_option = true;
	}
	else if (option == "save-ground-truth") {
		ground_truth_outfile_ = arg;
		is_valid_option = true;
	}
	else if (option == "log") {
		logfile_ = arg;
		std::ofstream logfile;
		logfile.open(logfile_);
		logfile.close();
		is_valid_option = true;
	}
	else if (option == "dnn-activation") {
		dnn_params_.activation_candidates.clear();
		std::stringstream activation_functions(arg);
		std::string activation;
		while (std::getline(activation_functions, activation, ',')) {
			if (activation == "SIGMOID") {
				dnn_params_.activation_candidates.push_back(FANN::activation_function_enum::SIGMOID);
			}
			else if (activation == "RELU") {
				dnn_params_.activation_candidates.push_back(FANN::activation_function_enum::RELU);
			}
			else {
				throw Error(std::string("Invalid argument ") + arg  + " for option dnn-activation. Usage: options = \"[--dnn-activation SIGMOID|RELU[,SIGMOID|RELU]] [...]\"");
			}
		}
		if (dnn_params_.activation_candidates.empty()) {
			throw Error(std::string("Invalid argument ") + arg  + " for option dnn-activation. Usage: options = \"[--dnn-activation SIGMOID|RELU[,SIGMOID|RELU]] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "dnn-hidden-layers-range") {
		std::stringstream hidden_layers_range(arg);
		std::string min_num_hidden_layers, max_num_hidden_layers;
		if (std::getline(hidden_layers_range, min_num_hidden_layers, ',') and std::getline(hidden_layers_range, max_num_hidden_layers, ',')) {
			try {
				dnn_params_.min_num_hidden_layers = std::stoi(min_num_hidden_layers);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option dnn-hidden-layers-range. Usage: options = \"[--dnn-hidden-layers-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
			try {
				dnn_params_.max_num_hidden_layers = std::stoi(max_num_hidden_layers);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option dnn-hidden-layers-range. Usage: options = \"[--dnn-hidden-layers-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
			if ((dnn_params_.min_num_hidden_layers > dnn_params_.max_num_hidden_layers) or (dnn_params_.min_num_hidden_layers < 0)) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option dnn-hidden-layers-range. Usage: options = \"[--dnn-hidden-layers-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option dnn-hidden-layers-range. Usage: options = \"[--dnn-hidden-layers-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "dnn-neurons-per-layer-range") {
		std::stringstream neurons_per_layer_range(arg);
		std::string min_neurons_per_layer, max_neurons_per_layer;
		if (std::getline(neurons_per_layer_range, min_neurons_per_layer, ',') and std::getline(neurons_per_layer_range, max_neurons_per_layer, ',')) {
			try {
				dnn_params_.min_num_neurons_per_layer = std::stoi(min_neurons_per_layer);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option dnn-neurons-per-layer-range. Usage: options = \"[--dnn-neurons-per-layer-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
			try {
				dnn_params_.max_num_neurons_per_layer = std::stoi(max_neurons_per_layer);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option dnn-neurons-per-layer-range. Usage: options = \"[--dnn-neurons-per-layer-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
			if ((dnn_params_.min_num_neurons_per_layer > dnn_params_.max_num_neurons_per_layer) or (dnn_params_.min_num_neurons_per_layer <= 0)) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option dnn-neurons-per-layer-range. Usage: options = \"[--dnn-neurons-per-layer-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option dnn-neurons-per-layer-range. Usage: options = \"[--dnn-neurons-per-layer-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "svm-gamma-exp-range") {
		std::stringstream gamma_exp_range(arg);
		std::string min_gamma_exp, max_gamma_exp;
		if (std::getline(gamma_exp_range, min_gamma_exp, ',') and std::getline(gamma_exp_range, max_gamma_exp, ',')) {
			try {
				svm_params_.min_gamma_exp = std::stoi(min_gamma_exp);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option svm-gamma-exp-range. Usage: options = \"[--svm-gamma-exp-range <smaller convertible to int>,<larger convertible to int>] [...]");
			}
			try {
				svm_params_.max_gamma_exp = std::stoi(max_gamma_exp);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option svm-gamma-exp-range. Usage: options = \"[--svm-gamma-exp-range <smaller convertible to int>,<larger convertible to int>] [...]");
			}
			if (svm_params_.min_gamma_exp > svm_params_.max_gamma_exp) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option svm-gamma-exp-range. Usage: options = \"[--svm-gamma-exp-range <smaller convertible to int>,<larger convertible to int>] [...]");
			}
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option svm-gamma-exp-range. Usage: options = \"[--svm-gamma-exp-range <smaller convertible to int>,<larger convertible to int>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "svm-c-exp-range") {
		std::stringstream c_exp_range(arg);
		std::string min_c_exp, max_c_exp;
		if (std::getline(c_exp_range, min_c_exp, ',') and std::getline(c_exp_range, max_c_exp, ',')) {
			try {
				svm_params_.min_c_exp = std::stoi(min_c_exp);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option svm-c-exp-range. Usage: options = \"[--svm-c-exp-range <smaller convertible to int>,<larger convertible to int>] [...]");
			}
			try {
				svm_params_.max_c_exp = std::stoi(max_c_exp);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option svm-c-exp-range. Usage: options = \"[--svm-c-exp-range <smaller convertible to int>,<larger convertible to int>] [...]");
			}
			if (svm_params_.min_c_exp > svm_params_.max_c_exp) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option svm-c-exp-range. Usage: options = \"[--svm-c-exp-range <smaller convertible to int>,<larger convertible to int>] [...]");
			}
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option svm-c-exp-range. Usage: options = \"[--svm-c-exp-range <smaller convertible to int>,<larger convertible to int>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "one-class-svm-likelihood") {
		if (arg == "TRUE") {
			one_class_svm_use_likelihood_ = true;
		}
		else if (arg == "FALSE") {
			one_class_svm_use_likelihood_ = false;
		}
		else {
			throw Error(std::string("Invalid argument ") + arg  + " for option one-class-svm-likelihood. Usage: options = \"[--one-class-svm-likelihood TRUE|FALSE] [...]\"");
		}
		is_valid_option = true;
	}
	if (dynamic_cast<IPFP<UserNodeLabel, UserEdgeLabel> *>(ground_truth_method_)) {
		if (ground_truth_options_ == "") {
			ground_truth_method_->set_options(std::string("--initial-solutions 80 --ratio-runs-from-initial-solutions 0.5 --lower-bound-method BRANCH_TIGHT --threads ") + std::to_string(this->num_threads_));
		}
		else {
			ground_truth_method_->set_options(ground_truth_options_ + " --threads " + std::to_string(this->num_threads_));
		}
	}
	else {
		if (ground_truth_options_ == "") {
			ground_truth_method_->set_options(std::string("--threads ") + std::to_string(this->num_threads_));
		}
		else {
			ground_truth_method_->set_options(ground_truth_options_ + " --threads " + std::to_string(this->num_threads_));
		}
	}
	is_valid_option = is_valid_option or ml_parse_option_(option, arg);
	return is_valid_option;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {
	ml_init_graph_(graph);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_default_post_graph_init_() {}

// === Definition of private helper functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
train_() {
	if (ml_method_ == DNN) {
		dnn_.train(dnn_training_data_, dnn_params_, outfile_, this->num_threads_);
	}
	else if (ml_method_ == SVM) {
		svm_.train(&svm_training_data_, svm_params_, num_features_, outfile_, this->num_threads_);
	}
	else {
		one_class_svm_.train(&svm_training_data_, one_class_svm_use_likelihood_, num_features_, outfile_);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
initialized_for_prediction_(const GEDGraph & g, const GEDGraph & h) const {
	return ((prediction_initialized_.first == g.id()) and (prediction_initialized_.second == h.id()));
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
generate_assignments_(const GEDGraph & g, const GEDGraph & h) {

	// Compute ground truth.
	NodeMap ground_truth(g.num_nodes(), h.num_nodes());
	if (compute_or_load_ground_truth_()) {
		if (ground_truth_infile_ != "") {
			this->ged_data_.load_node_map(ground_truth_infile_, g.id(), h.id(), ground_truth);
		}
		else {
			Result result;
			ground_truth_method_->run_as_util(g, h, result);
			ground_truth = result.node_map(0);
			if (ground_truth_outfile_ != "") {
				this->ged_data_.save_node_map(ground_truth_outfile_, g.id(), h.id(), ground_truth);
			}
		}
	}

	// Initialize variables used for computing the feature vectors
	ml_init_feature_variables_(g, h, this->num_threads_);
	std::vector<std::pair<std::size_t, std::size_t>> bad_assignments;
	std::size_t num_good_assignments{0};

	// Compute feature vectors. Skip bad assignments if called at initialization.
#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master <= g.num_nodes(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master <= h.num_nodes(); col_in_master++) {
			if ((row_in_master == g.num_nodes()) and (col_in_master == h.num_nodes())) {
				continue;
			}
			std::vector<double> feature_vector;
			bool good_assignment{false};
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				if (compute_or_load_ground_truth_()) {
					good_assignment = (ground_truth.image(row_in_master) == col_in_master);
				}
				if ((not this->initialized_) and (not good_assignment)) {
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						bad_assignments.emplace_back(row_in_master, col_in_master);
					}
					continue;
				}
				ml_populate_substitution_feature_vector_(g, h, row_in_master, col_in_master, feature_vector);
			}
			else if (row_in_master < g.num_nodes()) {
				if (compute_or_load_ground_truth_()) {
					good_assignment = (ground_truth.image(row_in_master) == GEDGraph::dummy_node());
				}
				if ((not this->initialized_) and (not good_assignment)) {
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						bad_assignments.emplace_back(row_in_master, col_in_master);
					}
					continue;
				}
				ml_populate_deletion_feature_vector_(g, row_in_master, feature_vector);
			}
			else {
				if (compute_or_load_ground_truth_()) {
					good_assignment = (ground_truth.pre_image(col_in_master) == GEDGraph::dummy_node());
				}
				if ((not this->initialized_) and (not good_assignment)) {
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						bad_assignments.emplace_back(row_in_master, col_in_master);
					}
					continue;
				}
				ml_populate_insertion_feature_vector_(h, col_in_master, feature_vector);
			}
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				assignments_.emplace_back(row_in_master, col_in_master, good_assignment, feature_vector);
				if (not this->initialized_) {
					num_good_assignments++;
				}
			}
		}
	}

	// Sample bad assignments if called at initialization and two class SVM or DNN is used.
	if ((not this->initialized_) and (ml_method_ != ONE_CLASS_SVM)) {
		std::random_device rng;
		std::mt19937 urng(rng());
		std::shuffle(bad_assignments.begin(), bad_assignments.end(), urng);
		std::size_t num_bad_assignments{bad_assignments.size()};
#ifdef _OPENMP
		omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
		for (std::size_t i = 0; i < std::min(num_good_assignments, num_bad_assignments); i++) {
			std::size_t row_in_master{bad_assignments.at(i).first};
			std::size_t col_in_master{bad_assignments.at(i).second};
			std::vector<double> feature_vector;
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				ml_populate_substitution_feature_vector_(g, h, row_in_master, col_in_master, feature_vector);
			}
			else if (row_in_master < g.num_nodes()) {
				ml_populate_deletion_feature_vector_(g, row_in_master, feature_vector);
			}
			else {
				ml_populate_insertion_feature_vector_(h, col_in_master, feature_vector);
			}
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				assignments_.emplace_back(row_in_master, col_in_master, false, feature_vector);
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
save_training_data_() {
	std::ofstream training_data_file(training_outfile_);
	for (Assignment_ assignment : assignments_) {
		training_data_file << assignment.to_string() << "\n";
	}
	training_data_file.close();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
load_or_generate_training_data_() {

	// Clear variables.
	assignments_.clear();
	dnn_feature_vectors_.clear();
	dnn_types_.clear();
	svm_feature_vectors_.clear();
	svm_types_.clear();


	// Load or generate training data.
	if (training_infile_ != "") {
		std::cout << "Loading training data from file " << training_infile_ << " ... " << std::flush;
		std::ifstream training_data_file(training_infile_);
		if (not training_data_file.good()) {
			throw Error("Error loading training data from file " + training_infile_ + ". File cannot be opened.");
		}
		std::string line;
		while(std::getline(training_data_file, line)) {
			assignments_.emplace_back(line, num_features_);
		}
		training_data_file.close();
		std::cout << "Done. Training data has size " << assignments_.size() << ".\n";
	}
	else {
		// Generate the assignments.
		ProgressBar progress_bar(this->ged_data_.num_graphs() * this->ged_data_.num_graphs());
		std::cout << "\rGenerating training data: " << progress_bar << std::flush;
		for (auto g = this->ged_data_.begin(); g != this->ged_data_.end(); g++) {
			if (this->ged_data_.is_shuffled_graph_copy(g->id())) {
				continue;
			}
			for (auto h = this->ged_data_.begin(); h != this->ged_data_.end(); h++) {
				if (this->ged_data_.is_shuffled_graph_copy(h->id())) {
					continue;
				}
				if (g->num_nodes() == 0 or h->num_nodes() == 0) {
					progress_bar.increment();
					std::cout << "\rGenerating training data: " << progress_bar << std::flush;
					continue;
				}
				// Generate assignments between g and h.
				if (this->ged_data_.shuffled_graph_copies_available() and (g->id() == h->id())) {
					generate_assignments_(*g, this->ged_data_.graph(this->ged_data_.id_shuffled_graph_copy(h->id())));
				}
				else {
					generate_assignments_(*g, *h);
				}
				progress_bar.increment();
				std::cout << "\rGenerating training data: " << progress_bar << std::flush;
			}
		}
		std::cout << "\n";
		std::random_shuffle(assignments_.begin(), assignments_.end());
	}

	// Initialize the specialized training data for DNN_ or SVM_.
	if (ml_method_ == DNN) {
		for (std::size_t pos{0}; pos < assignments_.size(); pos++) {
			dnn_feature_vectors_.push_back(assignments_[pos].dnn_feature_vector());
			dnn_types_.push_back(assignments_[pos].type());
		}
		dnn_training_data_.set_train_data(static_cast<unsigned int>(dnn_feature_vectors_.size()), static_cast<unsigned int>(num_features_), dnn_feature_vectors_.data(), 1, dnn_types_.data());
	}
	else {
		for (std::size_t pos{0}; pos < assignments_.size(); pos++) {
			svm_feature_vectors_.push_back(assignments_[pos].svm_feature_vector());
			if (ml_method_ == SVM) {
				svm_types_.push_back((*assignments_[pos].type()) + 1);
			}
			else {
				svm_types_.push_back(1);
			}
		}
		svm_training_data_.l = static_cast<int>(svm_feature_vectors_.size());
		svm_training_data_.y = svm_types_.data();
		svm_training_data_.x = svm_feature_vectors_.data();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
decision_value_(Assignment_ & assignment) {
	if (ml_method_ == DNN) {
		return dnn_.decision_value(assignment.dnn_feature_vector());
	}
	else if (ml_method_ == SVM) {
		return svm_.decision_value(assignment.svm_feature_vector());
	}
	return one_class_svm_.decision_value(assignment.svm_feature_vector());
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
load_config_file_() const {
	return (infile_ != "");
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
log_prediction_ratios_() const {
	return (logfile_ != "");
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
compute_or_load_ground_truth_() const {
	return (((not load_config_file_()) and (not this->initialized_)) or log_prediction_ratios_());
}

// === Definition of private class Assignment_. ===
template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
Assignment_(std::size_t row_in_master, std::size_t col_in_master, bool good_assignment, const std::vector<double> & feature_vector) :
row_in_master_{row_in_master},
col_in_master_{col_in_master},
type_{good_assignment ? 0.0 : 1.0},
dnn_feature_vector_(feature_vector),
svm_feature_vector_() {
	svm_node node;
	for (std::size_t index{0}; index < feature_vector.size(); index++) {
		if (feature_vector.at(index) != 0.0) {
			node.index = index + 1;
			node.value = feature_vector.at(index);
			svm_feature_vector_.push_back(node);
		}
	}
	node.index = -1;
	node.value = -1;
	svm_feature_vector_.push_back(node);
}

template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
Assignment_(const Assignment_ & assignment) :
row_in_master_{assignment.row_in_master_},
col_in_master_{assignment.col_in_master_},
type_{assignment.type_},
dnn_feature_vector_(assignment.dnn_feature_vector_),
svm_feature_vector_(assignment.svm_feature_vector_) {}

template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
Assignment_(const std::string & line, std::size_t num_features) :
row_in_master_{undefined()},
col_in_master_{undefined()},
type_{1.0},
dnn_feature_vector_(),
svm_feature_vector_() {
	std::istringstream line_stream(line);
	std::string type_str;
	if (not std::getline(line_stream, type_str, ' ')) {
		throw Error("Reading training data failed.\nExpected format of lines: \"0|1 <value 1> ... <value num_features>\"\nLine: \"" + line + "\"");
	}
	try {
		type_ = std::stod(type_str);
	}
	catch (...) {
		throw Error("Reading training data failed.\nExpected format of lines: \"0|1 <value 1> ... <value num_features>\"\nLine: \"" + line + "\"");
	}
	if ((type_ != 0) and (type_ != 1)) {
		throw Error("Reading training data failed.\nExpected format of lines: \"0|1 <value 1> ... <value num_features>\"\nLine: \"" + line + "\"");
	}
	std::string value_str;
	double value;
	int svm_index{0};
	svm_node node;
	while (std::getline(line_stream, value_str, ' ')) {
		svm_index++;
		try {
			value = std::stod(value_str);
		}
		catch (...) {
			throw Error("Reading training data failed.\nExpected format of lines: \"0|1 <value 1> ... <value num_features>\"\nLine: \"" + line + "\"");
		}
		dnn_feature_vector_.push_back(value);
		if (value != 0) {
			node.index = svm_index;
			node.value = value;
			svm_feature_vector_.push_back(node);
		}
	}
	node.index = -1;
	svm_feature_vector_.push_back(node);
	if (dnn_feature_vector_.size() != num_features) {
		throw Error("Reading training data failed.\nExpected format of lines: \"0|1 <value 1> ... <value num_features>\"\nLine: \"" + line + "\"");
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
row_in_master() const {
	return row_in_master_;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
col_in_master() const {
	return col_in_master_;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
to_string() const {
	std::string return_str(std::to_string(type_));
	for (double feature : dnn_feature_vector_) {
		return_str += " " + std::to_string(feature);
	}
	return return_str;
}

template<class UserNodeLabel, class UserEdgeLabel>
double *
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
dnn_feature_vector() {
	return dnn_feature_vector_.data();
}

template<class UserNodeLabel, class UserEdgeLabel>
struct svm_node *
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
svm_feature_vector() {
	return svm_feature_vector_.data();
}

template<class UserNodeLabel, class UserEdgeLabel>
double *
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
type() {
	return &type_;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
is_good_assignment() const {
	return (type_ == 0);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
Assignment_ ::
num_features() const {
	return dnn_feature_vector_.size();
}

// === Definition of private struct DNNParams_. ===
template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
DNNParams_ ::
DNNParams_() :
activation_candidates{FANN::activation_function_enum::RELU, FANN::activation_function_enum::SIGMOID},
min_num_hidden_layers{1},
max_num_hidden_layers{10},
min_num_neurons_per_layer{1},
max_num_neurons_per_layer{20} {}

// === Definition of private class DNN_. ===
template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
DNN_ ::
DNN_() :
neural_net_() {}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
DNN_ ::
load(const std::string & filename) {
	neural_net_.create_from_file(filename);
	return static_cast<std::size_t>(neural_net_.get_num_input());
}

template<class UserNodeLabel, class UserEdgeLabel>
float
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
DNN_ ::
cross_validate_(FANN::training_data & training_data, const MLBasedMethod::DNNParams_ & params, unsigned int num_hidden_layers, unsigned int num_neurons_per_layer, FANN::activation_function_enum hidden_activation) {

	std::size_t num_folds{5};
	std::size_t max_num_epochs{100};

	// Setup the neural network.
	std::vector<unsigned int> structure_dnn{training_data.num_input_train_data()};
	for (unsigned int hidden_layer{0}; hidden_layer < num_hidden_layers; hidden_layer++) {
		structure_dnn.push_back(num_neurons_per_layer);
	}
	structure_dnn.push_back(training_data.num_output_train_data());
	FANN::neural_net neural_net;
	neural_net.create_standard_array(structure_dnn.size(), structure_dnn.data());
	neural_net.set_activation_function_hidden(hidden_activation);
	neural_net.set_activation_steepness_hidden(1.0);
	neural_net.set_activation_function_output(FANN::activation_function_enum::SIGMOID);
	neural_net.set_activation_steepness_output(1.0);
	neural_net.set_train_error_function(FANN::error_function_enum::ERRORFUNC_LINEAR);
	neural_net.set_training_algorithm(FANN::training_algorithm_enum::TRAIN_INCREMENTAL);

	// Divide the training data into folds.
	unsigned int size_training_data{training_data.length_train_data()};
	unsigned int size_folds{size_training_data / static_cast<unsigned int>(num_folds)};
	double ** input_train_data{training_data.get_input()};
	double ** output_train_data{training_data.get_output()};
	std::vector<std::vector<double *>> folds_input(num_folds);
	std::vector<std::vector<double *>> folds_output(num_folds);
	std::size_t pos{0};
	for (unsigned int fold{0}; fold < num_folds - 1; fold++) {
		for (unsigned int counter{0}; counter < size_folds; counter++) {
			folds_input[fold].push_back(input_train_data[pos]);
			folds_output[fold].push_back(output_train_data[pos++]);
		}
	}
	for (; pos < size_training_data; pos++) {
		folds_input[num_folds - 1].push_back(input_train_data[pos]);
		folds_output[num_folds - 1].push_back(output_train_data[pos++]);
	}

	// Train and validate on the folds.
	float validation_error{0.0};
	for (unsigned int validation_fold{0}; validation_fold < num_folds; validation_fold++) {
		std::vector<double *> training_folds_input;
		std::vector<double *> training_folds_output;
		for (unsigned int fold{0}; fold < num_folds; fold++) {
			if (fold != validation_fold) {
				for (auto ptr : folds_input[fold]) {
					training_folds_input.push_back(ptr);
				}
				for (auto ptr : folds_output[fold]) {
					training_folds_output.push_back(ptr);
				}
			}
		}
		FANN::training_data training_folds_data;
		training_folds_data.set_train_data(training_folds_input.size(), training_data.num_input_train_data(), training_folds_input.data(), training_data.num_output_train_data(), training_folds_output.data());
		FANN::training_data validation_fold_data;
		validation_fold_data.set_train_data(folds_input[validation_fold].size(), training_data.num_input_train_data(), folds_input[validation_fold].data(), training_data.num_output_train_data(), folds_output[validation_fold].data());
		validation_error += train_and_validate_(neural_net, training_folds_data, validation_fold_data, max_num_epochs);
	}
	return validation_error;
}

template<class UserNodeLabel, class UserEdgeLabel>
float
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
DNN_ ::
train_and_validate_(FANN::neural_net & neural_net, FANN::training_data & training_data, FANN::training_data & validation_data, std::size_t max_num_epochs) {
	float min_validation_error{std::numeric_limits<float>::infinity()};
	float current_validation_error{std::numeric_limits<float>::infinity()};
	for (std::size_t epoch{0}; epoch < max_num_epochs; epoch++) {
		neural_net.train_epoch(training_data);
		current_validation_error = neural_net.test_data(validation_data);
		min_validation_error = std::min(min_validation_error, current_validation_error);
		if ((current_validation_error < 0.01) or (current_validation_error / min_validation_error > 1.05)) {
			break;
		}
	}
	return current_validation_error;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
DNN_ ::
train(FANN::training_data & training_data, const MLBasedMethod::DNNParams_ & params, const std::string & filename, std::size_t num_threads) {

	training_data.scale_input_train_data(0.0, 1.0);

	// Carry out cross-validation to determine the network structure.
	unsigned int optimal_num_hidden_layers{params.min_num_hidden_layers};
	unsigned int optimal_num_neurons_per_layer{params.min_num_neurons_per_layer};
	FANN::activation_function_enum optimal_hidden_activation{params.activation_candidates.at(0)};
	std::size_t num_activation_candidates{params.activation_candidates.size()};
	std::size_t num_hidden_layers_canditates{1 + params.max_num_hidden_layers - params.min_num_hidden_layers};
	std::size_t num_neurons_per_layer_canditates{1 + params.max_num_neurons_per_layer - params.min_num_neurons_per_layer};
	std::size_t progress_bar_size{1};
	if ((num_activation_candidates * num_hidden_layers_canditates * num_neurons_per_layer_canditates) > 1) {
		progress_bar_size = 1 + (num_activation_candidates * num_hidden_layers_canditates * num_neurons_per_layer_canditates);
	}

	ProgressBar progress_bar(progress_bar_size);
	std::cout << "\rTraining DNN: " << progress_bar << std::flush;
	if ((num_activation_candidates * num_hidden_layers_canditates * num_neurons_per_layer_canditates) > 1) {
		float optimal_validation_error{std::numeric_limits<float>::infinity()};
#ifdef _OPENMP
		omp_set_num_threads(num_threads - 1);
#pragma omp parallel for schedule(dynamic) if(num_threads > 1)
#endif
		for (unsigned int num_hidden_layers = params.min_num_hidden_layers; num_hidden_layers <= params.max_num_hidden_layers; num_hidden_layers++) {
			for (unsigned int num_neurons_per_layer = params.min_num_neurons_per_layer; num_neurons_per_layer <= params.max_num_neurons_per_layer; num_neurons_per_layer++) {
				for (FANN::activation_function_enum hidden_activation : params.activation_candidates) {
					float validation_error{cross_validate_(training_data, params, num_hidden_layers, num_neurons_per_layer, hidden_activation)};
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						if (validation_error < optimal_validation_error) {
							optimal_num_hidden_layers = num_hidden_layers;
							optimal_num_neurons_per_layer = num_neurons_per_layer;
							optimal_validation_error = validation_error;
							optimal_hidden_activation = hidden_activation;
						}
						progress_bar.increment();
						std::cout << "\rTraining DNN: " << progress_bar << std::flush;
					}
				}
			}
		}
	}

	// Setup the neural network.
	std::vector<unsigned int> structure_dnn{training_data.num_input_train_data()};
	for (unsigned int hidden_layer{0}; hidden_layer < optimal_num_hidden_layers; hidden_layer++) {
		structure_dnn.push_back(optimal_num_neurons_per_layer);
	}
	structure_dnn.push_back(training_data.num_output_train_data());
	neural_net_.create_standard_array(structure_dnn.size(), structure_dnn.data());
	neural_net_.set_activation_function_hidden(optimal_hidden_activation);
	neural_net_.set_activation_steepness_hidden(1.0);
	neural_net_.set_activation_function_output(FANN::activation_function_enum::SIGMOID);
	neural_net_.set_activation_steepness_output(1.0);
	neural_net_.set_train_error_function(FANN::error_function_enum::ERRORFUNC_LINEAR);
	neural_net_.set_training_algorithm(FANN::training_algorithm_enum::TRAIN_INCREMENTAL);

	// Divide the training data into training and validation data.
	double ** input_data{training_data.get_input()};
	double ** output_data{training_data.get_output()};
	unsigned int size_data{training_data.length_train_data()};
	std::vector<double *> train_input;
	std::vector<double *> train_output;
	std::vector<double *> valid_input;
	std::vector<double *> valid_output;
	for (unsigned int pos{0}; pos < size_data; pos++) {
		if (pos <= size_data / 5) {
			valid_input.push_back(input_data[pos]);
			valid_output.push_back(output_data[pos]);
		}
		else {
			train_input.push_back(input_data[pos]);
			train_output.push_back(output_data[pos]);
		}
	}

	// Train the neural network.
	std::size_t max_num_epochs_training{5000};
	FANN::training_data train_data;
	train_data.set_train_data(train_input.size(), training_data.num_input_train_data(), train_input.data(), training_data.num_output_train_data(), train_output.data());
	FANN::training_data valid_data;
	valid_data.set_train_data(valid_input.size(), training_data.num_input_train_data(), valid_input.data(), training_data.num_output_train_data(), valid_output.data());
	float valiation_error{train_and_validate_(neural_net_, train_data, valid_data, max_num_epochs_training)};
	progress_bar.increment();
	std::cout << "\rTraining DNN: " << progress_bar << std::flush;
	std::cout << "\nNetwork structure: " << training_data.num_input_train_data() << " x ";
	for (unsigned int layer{0}; layer < optimal_num_hidden_layers; layer++) {
		std::cout << optimal_num_neurons_per_layer << " x ";
	}
	std::cout << training_data.num_output_train_data() << ". Hidden layer activation: ";
	if (optimal_hidden_activation == FANN::activation_function_enum::SIGMOID) {
		std::cout << "Sigmoid";
	}
	else {
		std::cout << "ReLu";
	}
	std::cout << ". Validation error: " << valiation_error << "\n";

	// Save the trained neural network.
	if (filename != "") {
		neural_net_.save(filename);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
DNN_ ::
decision_value(double * feature_vector) {
	double * output = neural_net_.run(feature_vector);
	return (*output);
}

// === Definition of private struct SVMParams_. ===
template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
SVMParams_ ::
SVMParams_() :
min_gamma_exp{-3},
max_gamma_exp{3},
min_c_exp{-3},
max_c_exp{3},
min_nu{0.5},
max_nu{0.5} {}

// === Definition of private class SVM_. ===
template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
SVM_ ::
SVM_() :
svm_model_{nullptr} {}

template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
SVM_ ::
~SVM_() {
	svm_free_and_destroy_model(&svm_model_);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
SVM_ ::
load(const std::string & filename) {
	svm_model_ = svm_load_model(filename.c_str());
	if (svm_check_probability_model(svm_model_) == 0) {
		throw Error("SVM model does not support probability estimates.");
	}
	std::map<std::string, std::string> options;
	util::parse_config_file(filename + ".nf", options);
	return std::stoul(options.at("num_features"));
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
SVM_ ::
train(struct svm_problem * training_data, const MLBasedMethod::SVMParams_ & params, std::size_t num_features, const std::string & filename, std::size_t num_threads) {

	// Set the meta-parameters.
	struct svm_parameter svm_params;
	svm_params.gamma = std::pow(10, params.min_gamma_exp);
	svm_params.C = std::pow(10, params.min_c_exp);
	svm_params.svm_type = C_SVC;
	svm_params.kernel_type = RBF;
	svm_params.coef0 = 0;
	svm_params.degree = 0;
	svm_params.eps = 0.001;
	svm_params.cache_size = 100;
	svm_params.shrinking = 0;
	svm_params.probability = 1;
	svm_params.nr_weight = 0;
	svm_params.weight_label = nullptr;
	svm_params.weight = nullptr;
	svm_params.p = 0;
	svm_params.nu = 0;
	const char * error_msg;
	std::string error_msg_string;

	std::size_t progress_bar_size{1};
	if ((params.min_gamma_exp < params.max_gamma_exp) or (params.min_c_exp < params.max_c_exp)) {
		progress_bar_size += (1 + params.max_gamma_exp - params.min_gamma_exp) * (1 + params.max_c_exp - params.min_c_exp);
	}
	ProgressBar progress_bar(progress_bar_size);
	std::cout << "\rTraining SVM: " << progress_bar << std::flush;
	// Cross-validate to find the best values for gamma and C.
	if ((params.min_gamma_exp < params.max_gamma_exp) or (params.min_c_exp < params.max_c_exp)) {
		std::vector<double> predicted_y(training_data->l);
		std::size_t max_correct{0};
		int best_gamma_exp{params.min_gamma_exp};
		int best_c_exp{params.min_c_exp};
		std::size_t num_folds{5};
#ifdef _OPENMP
		omp_set_num_threads(num_threads - 1);
#pragma omp parallel for if(num_threads > 1)
#endif
		for (int gamma_exp = params.min_gamma_exp; gamma_exp <= params.max_gamma_exp; gamma_exp++) {
			struct svm_parameter local_svm_params = svm_params;
			local_svm_params.probability = 0;
			local_svm_params.gamma = std::pow(10, static_cast<double>(gamma_exp));
			for (int c_exp{params.min_c_exp}; c_exp <= params.max_c_exp; c_exp++) {
				local_svm_params.C = std::pow(10, static_cast<double>(c_exp));
				error_msg = svm_check_parameter(training_data, &local_svm_params);
				if (error_msg) {
					error_msg_string = std::string(error_msg);
					throw Error(error_msg_string);
				}
				svm_cross_validation(training_data, &local_svm_params, num_folds, predicted_y.data());
				std::size_t correct{0};
				for (int i{0}; i < training_data->l; i++) {
					if (training_data->y[i] == predicted_y.at(i)) {
						correct++;
					}
				}
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					if (correct > max_correct) {
						max_correct = correct;
						best_gamma_exp = gamma_exp;
						best_c_exp = c_exp;
					}
					progress_bar.increment();
					std::cout << "\rTraining SVM: " << progress_bar << std::flush;
				}
			}
		}
		svm_params.gamma = std::pow(10, static_cast<double>(best_gamma_exp));
		svm_params.C = std::pow(10, static_cast<double>(best_c_exp));
	}

	// Train the SVM with the found parameters.
	error_msg = svm_check_parameter(training_data, &svm_params);
	if (error_msg) {
		error_msg_string = std::string(error_msg);
		throw Error(error_msg_string);
	}
	svm_model_ = svm_train(training_data, &svm_params);
	progress_bar.increment();
	std::cout << "\rTraining SVM: " << progress_bar << std::flush << "\n";
	if (svm_check_probability_model(svm_model_) == 0) {
		throw Error("SVM model does not support probability estimates.");
	}

	// Save the trained SVM.
	if (filename != "") {
		if (svm_save_model(filename.c_str(), svm_model_)) {
			throw Error("Cannot save SVM model to " + filename + ".");
		}
		std::map<std::string, std::string> options;
		options["num_features"] = std::to_string(num_features);
		util::save_as_config_file(filename + ".nf", options);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
SVM_ ::
decision_value(struct svm_node * feature_vector) const {
	std::vector<double> probability_estimates(2);
	double predicted_type{svm_predict_probability(svm_model_, feature_vector, probability_estimates.data())};
	double higher_prob = std::max(probability_estimates.at(0), probability_estimates.at(1));
	double lower_prob = std::min(probability_estimates.at(0), probability_estimates.at(1));
	if (predicted_type == 2) {
		return higher_prob;
	}
	else if (predicted_type != 1) {
		throw Error("Unexpected predicted type " + std::to_string(predicted_type) + ".");
	}
	return lower_prob;
}

// === Definition of private class OneClassSVM_. ===
template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
OneClassSVM_ ::
OneClassSVM_() :
svm_model_{nullptr},
rho_{0},
sum_alpha_{0},
scale_factor_{1},
use_likelihood_{false} {}

template<class UserNodeLabel, class UserEdgeLabel>
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
OneClassSVM_ ::
~OneClassSVM_() {
	svm_free_and_destroy_model(&svm_model_);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
OneClassSVM_ ::
load(const std::string & filename, bool use_likelihood) {
	use_likelihood_ = use_likelihood;
	svm_model_ = svm_load_model(filename.c_str());
	std::map<std::string, std::string> options;
	util::parse_config_file(filename + ".nf", options);
	std::size_t num_features{std::stoul(options.at("num_features"))};
	compute_rho_and_scale_factor_(num_features);
	return num_features;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
OneClassSVM_ ::
train(struct svm_problem * training_data, bool use_likelihood, std::size_t num_features, const std::string & filename) {

	// Set the meta-parameters.
	use_likelihood_ = use_likelihood;
	struct svm_parameter svm_params;
	svm_params.gamma = 1.0 / static_cast<double>(num_features);
	svm_params.C = 0;
	svm_params.svm_type = ONE_CLASS;
	svm_params.kernel_type = RBF;
	svm_params.coef0 = 0;
	svm_params.degree = 0;
	svm_params.eps = 0.001;
	svm_params.cache_size = 100;
	svm_params.shrinking = 0;
	svm_params.probability = 0;
	svm_params.nr_weight = 0;
	svm_params.weight_label = nullptr;
	svm_params.weight = nullptr;
	svm_params.p = 0;
	svm_params.nu = 0.5;
	const char * error_msg;
	std::string error_msg_string;

	std::size_t progress_bar_size{1};
	ProgressBar progress_bar(progress_bar_size);
	std::cout << "\rTraining one class SVM: " << progress_bar << std::flush;

	// Train the SVM.
	error_msg = svm_check_parameter(training_data, &svm_params);
	if (error_msg) {
		error_msg_string = std::string(error_msg);
		throw Error(error_msg_string);
	}
	svm_model_ = svm_train(training_data, &svm_params);
	progress_bar.increment();
	std::cout << "\rTraining SVM: " << progress_bar << std::flush << "\n";
	compute_rho_and_scale_factor_(num_features);

	// Save the trained SVM.
	if (filename != "") {
		if (svm_save_model(filename.c_str(), svm_model_)) {
			throw Error("Cannot save SVM model to " + filename + ".");
		}
		std::map<std::string, std::string> options;
		options["num_features"] = std::to_string(num_features);
		util::save_as_config_file(filename + ".nf", options);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
OneClassSVM_ ::
decision_value(struct svm_node * feature_vector) const {
	std::vector<double> dec_values(1);
	double predicted_type{svm_predict_values(svm_model_, feature_vector, dec_values.data())};
	if (predicted_type != 1 and predicted_type != -1) {
		throw Error("Unexpected predicted type " + std::to_string(predicted_type) + ".");
	}
	if (use_likelihood_) {
		return (1 - ((dec_values.at(0) - rho_) * scale_factor_));
	}
	else if (dec_values.at(0) >= 0.0) {
		return (0.5 - 0.5 * (dec_values.at(0) / (sum_alpha_ - rho_)));
	}
	return (0.5 - 0.5 * (dec_values.at(0) / rho_));
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
OneClassSVM_ ::
compute_rho_and_scale_factor_(std::size_t num_features) {
	rho_ = svm_model_->rho[0];
	double * alpha{svm_model_->sv_coef[0]};
	sum_alpha_ = 0;
	for (int i{0}; i < svm_model_->l; i++) {
		sum_alpha_ += alpha[i];
	}
	scale_factor_ = std::pow(svm_model_->param.gamma / pi(), static_cast<double>(num_features) / 2.0) / sum_alpha_;
}

// === Default definitions of virtual member functions to be overridden by derived classes. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_init_feature_variables_(const GEDGraph & g, const GEDGraph & h, std::size_t num_threads) {}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_valid_options_string_() const {
	return "";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_parse_option_(const std::string & option, const std::string & arg) {
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_set_default_options_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_init_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_init_graph_(const GEDGraph & graph) {}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_get_num_features_() {
	return undefined();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_init_for_num_features_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_populate_substitution_feature_vector_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k, std::vector<double> & feature_vector) {}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_populate_deletion_feature_vector_(const GEDGraph & g, GEDGraph::NodeID i, std::vector<double> & feature_vector) {}

template<class UserNodeLabel, class UserEdgeLabel>
void
MLBasedMethod<UserNodeLabel, UserEdgeLabel>::
ml_populate_insertion_feature_vector_(const GEDGraph & h, GEDGraph::NodeID k, std::vector<double> & feature_vector) {}

}

#endif /* SRC_METHODS_ML_BASED_METHOD_IPP_ */
