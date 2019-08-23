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
 * @file  ml_based_method.hpp
 * @brief ged::MLBasedMethod class declaration.
 */

#ifndef SRC_METHODS_ML_BASED_METHOD_HPP_
#define SRC_METHODS_ML_BASED_METHOD_HPP_

namespace ged {

/*!
 * @brief Abstract class for methods that transform GED to LSAPE by using a SVM or a DNN to predict the cost of node edit operations.
 * @details All derived classes support the following options in addition to the ones supported by ged::LSAPEBasedMethod:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--load @<filename@></tt> | path to existing configuration file | not specified | must be specified, if the method is run without prior initialization |
 * | <tt>\--save @<filename@></tt> | path where to save configuration file | not specified | n.a. |
 * | <tt>\--load-train @<filename@></tt> | path to existing training data file | not specified | n.a. |
 * | <tt>\--save-train @<filename@></tt> | path where to save training data | not specified | n.a. |
 * | <tt>\--load-ground-truth @<filename@></tt> | path to existing ground truth data file | not specified | n.a. |
 * | <tt>\--save-ground-truth @<filename@></tt> | path where to save ground truth | not specified | n.a. |
 * | <tt>\--log @<filename@></tt> | path where to save log data | not specified | specify this option for output of accuracy, precision, and recall |
 * | <tt>\--ml-method DNN\|SVM\|ONE_CLASS_SVM</tt> | employed machine learning method | @p DNN | if @p DNN, a fully connected neural network with one output neuron is used and the options @p \--svm-@<suffix@> and @p \--one-class-svm-likelihood have no effect <br> if @p SVM, a two-class support vector machine with RBF kernel is used and the options @p \--dnn-@<suffix@> and @p \--one-class-svm-likelihood have no effect <br> if @p ONE_CLASS_SVM, a one class support vector machine with RBF kernel is used and the options @p \--dnn-@<suffix@> and @p \--svm-@<suffix@> have no effect |
 * | <tt>\--ground-truth-method ANCHOR_AWARE_GED\|%F1\|%F2\|COMPACT_MIP\|%IPFP</tt> | method for computing the ground truth | @p %IPFP | the methods %F1, %F2, and COMPACT_MIP are available only if GEDLIB is installed with Gurobi |
 * | <tt>\--ground-truth-options '[--@<option@> @<arg@>] [...]'</tt> | options string passed to the ground truth method | @p '' | ged::AnchorAwareGED, ged::F1, ged::F2, ged::CompactMIP, ged::IPFP |
 * | <tt>\--dnn-activation SIGMOID\|RELU[,SIGMOID\|RELU]</tt> | activation functions for hidden layers | <tt>SIGMOID,RELU</tt> | if more than one function is specified, the best choice is determined via cross-validation during training |
 * | <tt>\--dnn-hidden-layers-range @<smaller convertible to int greater 0@>,@<larger convertible to int greater 0@></tt> | range that specifies possible number of hidden layers | <tt>1,10</tt> | if the range is larger than one, the best choice is determined via cross-validation during training |
 * | <tt>\--dnn-neurons-per-layer-range @<smaller convertible to int greater 0@>,@<larger convertible to int greater 0@></tt> | range that specifies possible number of neurons per hidden layer | <tt>1,20</tt> | if the range is larger than one, the best choice is determined via cross-validation during training |
 * | <tt>\--svm-gamma-exp-range @<smaller convertible to int@>,@<larger convertible to int@></tt> | range that specifies possible exponents to the basis 10 of the parameter @f$\gamma@f$ used by SVM | <tt>-3,4</tt> | if the range is larger than one, the best choice is determined via cross-validation during training |
 * | <tt>\--svm-c-exp-range @<smaller convertible to int@>,@<larger convertible to int@></tt> | range that specifies possible exponents to the basis 10 of the parameter @f$C@f$ used by SVM | <tt>-3,4</tt> | if the range is larger than one, the best choice is determined via cross-validation during training |
 * | <tt>\--one-class-svm-likelihood TRUE\|FALSE</tt> | use likelihood to define probability estimates based on one class SVM output | @p TRUE | n.a. |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class MLBasedMethod : public LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~MLBasedMethod() = 0;

	MLBasedMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

	/*!
	 * @brief Predicts the type of a node assignment.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[in] assignment An assignment of a node in @p g or ged::dummy_node() to a node in @p h or ged::dummy_node().
	 * @return A value between 0 and 1. If smaller than 0.5, the assignment is predicted to be good, otherwise it is predicted to be bad.
	 */
	double predict(const GEDGraph & g, const GEDGraph & h, const NodeMap::Assignment & assignment);


protected:

	/*!
	 * @brief The size of the feature vectors.
	 */
	std::size_t num_features_;

private:

	enum MLMethod_ {DNN, SVM, ONE_CLASS_SVM};

	class Assignment_ {

	public:

		Assignment_(std::size_t row_in_master, std::size_t col_in_master, bool good_assignment, const std::vector<double> & feature_vector);

		Assignment_(const std::string & line, std::size_t num_features);

		Assignment_(const Assignment_ & assignment);

		std::size_t row_in_master() const;

		std::size_t col_in_master() const;

		std::string to_string() const;

		double * dnn_feature_vector();

		struct svm_node * svm_feature_vector();

		double * type();

		bool is_good_assignment() const;

		std::size_t num_features() const;

	private:

		std::size_t row_in_master_;

		std::size_t col_in_master_;

		double type_;

		std::vector<double> dnn_feature_vector_;

		std::vector<struct svm_node> svm_feature_vector_;
	};

	struct DNNParams_ {

		DNNParams_();

		std::vector<FANN::activation_function_enum> activation_candidates;

		unsigned int min_num_hidden_layers;

		unsigned int max_num_hidden_layers;

		unsigned int min_num_neurons_per_layer;

		unsigned int max_num_neurons_per_layer;
	};

	class DNN_ {

	public:

		DNN_();

		std::size_t load(const std::string & filename);

		void train(FANN::training_data & training_data, const MLBasedMethod::DNNParams_ & params, const std::string & filename, std::size_t num_threads);

		double decision_value(double * feature_vector);

	private:

		float cross_validate_(FANN::training_data & training_data, const MLBasedMethod::DNNParams_ & params, unsigned int num_hidden_layers, unsigned int num_neurons_per_layer, FANN::activation_function_enum hidden_activation);

		float train_and_validate_(FANN::neural_net & neural_net, FANN::training_data & training_data, FANN::training_data & validation_data, std::size_t max_num_epochs);

		FANN::neural_net neural_net_;
	};

	struct SVMParams_ {

		SVMParams_();

		int min_gamma_exp;

		int max_gamma_exp;

		int min_c_exp;

		int max_c_exp;

		double min_nu;

		double max_nu;
	};

	class SVM_ {

	public:

		~SVM_();

		SVM_();

		std::size_t load(const std::string & filename);

		void train(struct svm_problem * training_data, const MLBasedMethod::SVMParams_ & params, std::size_t num_features, const std::string & filename, std::size_t num_threads);

		double decision_value(struct svm_node * feature_vector) const;

	private:

		struct svm_model * svm_model_;

	};

	class OneClassSVM_ {

	public:

		~OneClassSVM_();

		OneClassSVM_();

		std::size_t load(const std::string & filename, bool use_likelihood);

		void train(struct svm_problem * training_data, bool use_likelihood, std::size_t num_features, const std::string & filename);

		double decision_value(struct svm_node * feature_vector) const;

	private:

		struct svm_model * svm_model_;

		double rho_;

		double sum_alpha_;

		double scale_factor_;

		bool use_likelihood_;

		void compute_rho_and_scale_factor_(std::size_t num_features);
	};

	std::pair<GEDGraph::GraphID, GEDGraph::GraphID> prediction_initialized_;

	std::vector<Assignment_> assignments_;

	MLMethod_ ml_method_;

	GEDMethod<UserNodeLabel, UserEdgeLabel> * ground_truth_method_;

	std::string ground_truth_options_;

	DNNParams_ dnn_params_;

	DNN_ dnn_;

	FANN::training_data dnn_training_data_;

	std::vector<double *> dnn_feature_vectors_;

	std::vector<double *> dnn_types_;

	SVMParams_ svm_params_;

	SVM_ svm_;

	struct svm_problem svm_training_data_;

	std::vector<struct svm_node *> svm_feature_vectors_;

	std::vector<double> svm_types_;

	bool one_class_svm_use_likelihood_;

	OneClassSVM_ one_class_svm_;

	std::string infile_;

	std::string outfile_;

	std::string logfile_;

	std::string training_infile_;

	std::string training_outfile_;

	std::string ground_truth_infile_;

	std::string ground_truth_outfile_;

	// Member functions inherited from LSAPEBasedMethod.

	virtual void lsape_init_() final;

	virtual void lsape_pre_graph_init_(bool called_at_runtime) final;

	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) final;

	virtual std::string lsape_valid_options_string_() const final;

	virtual void lsape_set_default_options_() final;

	virtual bool lsape_parse_option_(const std::string & option, const std::string & value) final;

	virtual void lsape_init_graph_(const GEDGraph & graph) final;

	virtual void lsape_default_post_graph_init_() final;

	// Private helper functions.

	bool initialized_for_prediction_(const GEDGraph & g, const GEDGraph & h) const;

	void generate_assignments_(const GEDGraph & g, const GEDGraph & h);

	void load_or_generate_training_data_();

	void save_training_data_();

	bool load_config_file_() const;

	bool log_prediction_ratios_() const;

	bool compute_or_load_ground_truth_() const;

	double decision_value_(Assignment_ & assignmment);

	void train_();

	// Virtual functions to be overridden by derived classes.

	/*!
	 * @brief Initializes variables that are used for populating the feature vectors of assignments between two input graphs.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[in] num_threads The number of available threads.
	 * @note Must be overridden by derived classes that require initialization of variables used for populating the feature vectors of assignments between two input graphs.
	 */
	virtual void ml_init_feature_variables_(const GEDGraph & g, const GEDGraph & h, std::size_t num_threads);

	/*!
	 * @brief Returns string of all valid options that are not among the ones shared by all derived classes of ged::MLBasedMethod.
	 * @return String of the form @"[--@<option@> @<arg@>] [...]@".
	 * @note Must be overridden by derived classes that have options that are not among the ones shared by all derived classes of ged::MLBasedMethod.
	 */
	virtual std::string ml_valid_options_string_() const;

	/*!
	 * @brief Parses one option that is not among the ones shared by all derived classes of ged::MLBasedMethod.
	 * @param[in] option The name of the option.
	 * @param[in] arg The argument of the option.
	 * @return Returns true if @p option is a valid option name for the method and false otherwise.
	 * @note Must be overridden by derived classes that have options that are not among the ones shared by all derived classes of ged::MLBasedMethod.
	 */
	virtual bool ml_parse_option_(const std::string & option, const std::string & arg);

	/*!
	 * @brief Sets all options that are not among the ones shared by all derived classes of ged::MLBasedMethod to default values.
	 * @note Must be overridden by derived classes that have options that are not among the ones shared by all derived classes of ged::MLBasedMethod.
	 */
	virtual void ml_set_default_options_();

	/*!
	 * @brief Initializes the method after initializing the global variables for the graphs.
	 * @note Must be overridden by derived classes of ged::MLBasedMethod that require custom initialization.
	 */
	virtual void ml_init_();

	/*!
	 * @brief Initializes global variables for one graph.
	 * @param[in] graph Graph for which the global variables have to be initialized.
	 * @note Must be overridden by derived classes that require to initialize custom global variables.
	 */
	virtual void ml_init_graph_(const GEDGraph & graph);

	/*!
	 * @brief Returns the number of features.
	 * @return Size of the feature vectors employed by the derived class.
	 * @note Must be overridden by derived classes.
	 */
	virtual std::size_t ml_get_num_features_();

	/*!
	 * @brief Initializes the derived class for running with feature vectors of size ged::MLBasedMethod::num_features_.
	 * @note Must be overridden by derived classes that can be set up with features vectors of various sizes.
	 */
	virtual void ml_init_for_num_features_();

	/*!
	 * @brief Computes substitution feature vector.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[in] i ID of node in @p g that has to be substituted.
	 * @param[in] k ID of node in @p h that has to be substituted.
	 * @param[out] feature_vector Return variable. Empty when called. Must be of size ged::MLBasedMethod::num_features_ when exiting.
	 * @note Must be overridden by derived classes.
	 */
	virtual void ml_populate_substitution_feature_vector_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k, std::vector<double> & feature_vector);

	/*!
	 * @brief Computes deletion feature vector.
	 * @param[in] g Input graph.
	 * @param[in] i ID of node in @p g that has to be deleted.
	 * @param[out] feature_vector Return variable. Empty when called. Must be of size ged::MLBasedMethod::num_features_ when exiting.
	 * @note Must be overridden by derived classes.
	 */
	virtual void ml_populate_deletion_feature_vector_(const GEDGraph & g, GEDGraph::NodeID i, std::vector<double> & feature_vector);

	/*!
	 * @brief Computes insertion feature vector.
	 * @param[in] h Input graph.
	 * @param[in] k ID of node in @p h that has to be inserted.
	 * @param[out] feature_vector Return variable. Empty when called. Must be of size ged::MLBasedMethod::num_features_ when exiting.
	 * @note Must be overridden by derived classes.
	 */
	virtual void ml_populate_insertion_feature_vector_(const GEDGraph & h, GEDGraph::NodeID k, std::vector<double> & feature_vector);

};

}

#endif /* SRC_METHODS_ML_BASED_METHOD_HPP_ */

