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
 * @file bipartite_ml.ipp
 * @brief ged::BipartiteML class definition.
 */

#ifndef SRC_METHODS_BIPARTITE_ML_IPP_
#define SRC_METHODS_BIPARTITE_ML_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
BipartiteML<UserNodeLabel, UserEdgeLabel>::
~BipartiteML() {
	delete lsape_method_;
}

template<class UserNodeLabel, class UserEdgeLabel>
BipartiteML<UserNodeLabel, UserEdgeLabel>::
BipartiteML(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
MLBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
lsape_method_{new Bipartite<UserNodeLabel, UserEdgeLabel>(this->ged_data_)},
lsape_method_options_(""),
lsape_instance_(),
global_features_(),
row_features_(),
col_features_() {}

// === Definitions of member functions inherited from MLBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_init_() {
	lsape_method_->init();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_set_default_options_() {
	delete lsape_method_;
	lsape_method_ = new Bipartite<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
	lsape_method_options_ = std::string("");
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_valid_options_string_() const {
	return "[--lsape-method <arg>] [--lsape-options <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_parse_option_(const std::string & option, const std::string & arg) {
	bool is_valid_option{false};
	if (option == "lsape-method") {
		if (arg == "BRANCH_FAST") {
			lsape_method_ = new BranchFast<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BRANCH_UNIFORM") {
			lsape_method_ = new BranchUniform<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BRANCH") {
			lsape_method_ = new Branch<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "NODE") {
			lsape_method_ = new Node<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "RING") {
			lsape_method_ = new Ring<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "SUBGRAPH") {
			lsape_method_ = new Subgraph<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "WALKS") {
			lsape_method_ = new Walks<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg != "BIPARTITE") {
			throw Error("Invalid argument \"" + arg + "\" for option lsape-method. Usage: options = \"[--lsape-method BIPARTITE|BRANCH_FAST|BRANCH_UNIFORM|BRANCH|NODE|RING|SUBGRAPH|WALKS] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "lsape-options") {
		lsape_method_options_ = arg;
		std::size_t bad_option_start{lsape_method_options_.find("--threads")};
		std::size_t next_option_start;
		if (bad_option_start != std::string::npos) {
			next_option_start = lsape_method_options_.find("--", bad_option_start + 1);
			if (next_option_start != std::string::npos) {
				lsape_method_options_ = lsape_method_options_.substr(0, bad_option_start) + lsape_method_options_.substr(next_option_start);
			}
			else {
				lsape_method_options_ = lsape_method_options_.substr(0, bad_option_start);
			}
		}
		is_valid_option = true;
	}
	if (lsape_method_options_ != "") {
		lsape_method_->set_options(lsape_method_options_ + " --threads " + std::to_string(this->num_threads_));
	}
	else {
		lsape_method_->set_options(std::string("--threads ") + std::to_string(this->num_threads_));
	}
	return is_valid_option;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_init_feature_variables_(const GEDGraph & g, const GEDGraph & h, std::size_t num_threads) {
	populate_lsape_instance_(g, h, num_threads);
	Eigen::ArrayXXd substitution_matrix(lsape_instance_.topLeftCorner(g.num_nodes(), h.num_nodes()));
	compute_global_features_(substitution_matrix);
	row_features_.init(substitution_matrix);
	col_features_.init(substitution_matrix);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_populate_substitution_feature_vector_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k, std::vector<double> & feature_vector) {
	std::size_t row{i};
	std::size_t col{k};
	feature_vector.clear();
	add_global_features_(feature_vector);
	add_cell_features_(row, col, this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k)), feature_vector);
	row_features_.add_features_(lsape_instance_, row, col, feature_vector);
	col_features_.add_features_(lsape_instance_, row, col, feature_vector);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_populate_deletion_feature_vector_(const GEDGraph & g, GEDGraph::NodeID i, std::vector<double> & feature_vector) {
	std::size_t row{i};
	std::size_t col{static_cast<std::size_t>(lsape_instance_.cols() - 1)};
	feature_vector.clear();
	add_global_features_(feature_vector);
	add_cell_features_(row, col, this->ged_data_.node_cost(g.get_node_label(i), dummy_label()), feature_vector);
	row_features_.add_features_(lsape_instance_, row, col, feature_vector);
	col_features_.add_features_(lsape_instance_, row, col, feature_vector);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_populate_insertion_feature_vector_(const GEDGraph & h, GEDGraph::NodeID k, std::vector<double> & feature_vector) {
	std::size_t row{static_cast<std::size_t>(lsape_instance_.rows() - 1)};
	std::size_t col{k};
	feature_vector.clear();
	add_global_features_(feature_vector);
	add_cell_features_(row, col, this->ged_data_.node_cost(dummy_label(), h.get_node_label(k)), feature_vector);
	row_features_.add_features_(lsape_instance_, row, col, feature_vector);
	col_features_.add_features_(lsape_instance_, row, col, feature_vector);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ml_get_num_features_() {
	return 24;
}

// === Definition of private class RowFeatures_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BipartiteML<UserNodeLabel, UserEdgeLabel>::
RowFeatures_::
RowFeatures_() :
maxima_(),
minima_(),
means_(),
deviations_(),
leaders_(),
intervals_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
RowFeatures_::
init(const Eigen::ArrayXXd & substitution_matrix) {

	// Compute row maxima.
	maxima_ = substitution_matrix.rowwise().maxCoeff();

	// Compute row minima and their positions.
	minima_.resize(substitution_matrix.rows());
	std::vector<RowVector_::Index> col_min(substitution_matrix.rows());
	for (auto row = 0; row < substitution_matrix.rows(); row++) {
		minima_(row) = substitution_matrix.row(row).minCoeff(&col_min[row]);
	}

	// Compute row means.
	means_ = substitution_matrix.rowwise().mean();

	// Compute row deviations.
	if (substitution_matrix.rows() <= 1) {
		deviations_.resize(substitution_matrix.rows());
		for (auto row = 0; row < substitution_matrix.rows(); row++) {
			deviations_(row) = 0.0;
		}
	}
	else {
		deviations_ = ((substitution_matrix.colwise() - substitution_matrix.rowwise().mean()).square().rowwise().sum() / (substitution_matrix.rows() - 1)).sqrt();
	}

	// Compute row leaders.
	RowVector_ second_minima(substitution_matrix.rows());
	if (substitution_matrix.rows() == 1) {
		second_minima = minima_;
	}
	else {
		for (auto row = 0; row < substitution_matrix.rows(); row++) {
			second_minima(row) = maxima_(row);
			for (auto col = 0; col < substitution_matrix.cols(); col++) {
				if ((col != col_min.at(row)) and (substitution_matrix(row, col) < second_minima(row))) {
					second_minima(row) = substitution_matrix(row, col);
				}
			}
		}
	}
	leaders_ = (minima_ - second_minima);
	for (auto row = 0; row < substitution_matrix.rows(); row++) {
		if (second_minima(row) != 0) {
			leaders_(row) /= second_minima(row);
		}
	}

	// Compute row intervals.
	intervals_ = maxima_ - minima_;
	if (intervals_.mean() > 0) {
		intervals_ /= intervals_.mean();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
RowFeatures_::
add_features_(const Eigen::ArrayXXd & matrix, std::size_t row, std::size_t col, std::vector<double> & feature_vector) const {

	// Add zeroes if row is the dummy row.
	if (row == static_cast<std::size_t>(matrix.rows() - 1)) {
		for (unsigned short counter{0}; counter < 9; counter++) {
			feature_vector.push_back(0.0);
		}
		return;
	}

	// Add row maximum feature.
	feature_vector.push_back(maxima_(row));

	// Add row minimum feature.
	feature_vector.push_back(minima_(row));

	// Add row mean feature.
	feature_vector.push_back(means_(row));

	// Add row deviation feature.
	feature_vector.push_back(deviations_(row));

	// Add row uniqueness feature.
	double uniqueness{(matrix.row(row) - matrix(row, col)).abs().maxCoeff()};
	feature_vector.push_back(uniqueness);

	// Add row divergence feature.
	double divergence{(matrix.row(row) - matrix(row, col)).abs().sum()};
	if (((matrix.cols() - 1) * (matrix(row, col) - means_(row))) != 0) {
		divergence /= ((matrix.cols() - 1) * (matrix(row, col) - means_(row)));
	}
	feature_vector.push_back(divergence);

	// Add row leader feature.
	feature_vector.push_back(leaders_(row));

	// Add row interval feature.
	feature_vector.push_back(intervals_(row));

	// Add row outlierness feature.
	double outlierness{matrix(row, col)};
	if (means_(row) != deviations_(row)) {
		outlierness /= (means_(row) - deviations_(row));
	}
	feature_vector.push_back(outlierness);
}

// === Definition of private class ColFeatures_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ColFeatures_::
ColFeatures_() :
maxima_(),
minima_(),
means_(),
deviations_(),
leaders_(),
intervals_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ColFeatures_::
init(const Eigen::ArrayXXd & substitution_matrix) {

	// Compute column maxima.
	maxima_ = substitution_matrix.colwise().maxCoeff();

	// Compute column minima and their positions.
	minima_.resize(substitution_matrix.cols());
	std::vector<ColumnVector_::Index> row_min(substitution_matrix.cols());
	for (auto col = 0; col < substitution_matrix.cols(); col++) {
		minima_(col) = substitution_matrix.col(col).minCoeff(&row_min[col]);
	}

	// Compute column means.
	means_ = substitution_matrix.colwise().mean();

	// Compute column deviations.
	if (substitution_matrix.cols() <= 1) {
		deviations_.resize(substitution_matrix.cols());
		for (auto col = 0; col < substitution_matrix.cols(); col++) {
			deviations_(col) = 0.0;
		}
	}
	else {
		deviations_ = ((substitution_matrix.rowwise() - substitution_matrix.colwise().mean()).square().colwise().sum() / (substitution_matrix.cols() - 1)).sqrt();
	}

	// Compute column leaders.
	ColumnVector_ second_minima(substitution_matrix.cols());
	if (substitution_matrix.cols() == 1) {
		second_minima = minima_;
	}
	else {
		for (auto col = 0; col < substitution_matrix.cols(); col++) {
			second_minima(col) = maxima_(col);
			for (auto row = 0; row < substitution_matrix.rows(); row++) {
				if ((row != row_min.at(col)) and (substitution_matrix(row, col) < second_minima(col))) {
					second_minima(col) = substitution_matrix(row, col);
				}
			}
		}
	}
	leaders_ = (minima_ - second_minima);
	for (auto col = 0; col < substitution_matrix.cols(); col++) {
		if (second_minima(col) != 0) {
			leaders_(col) /= second_minima(col);
		}
	}

	// Compute row intervals.
	intervals_ = maxima_ - minima_;
	if (intervals_.mean() > 0) {
		intervals_ /= intervals_.mean();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
ColFeatures_::
add_features_(const Eigen::ArrayXXd & matrix, std::size_t row, std::size_t col, std::vector<double> & feature_vector) const {

	// Add zeroes if col is the dummy column.
	if (col == static_cast<std::size_t>(matrix.cols() - 1)) {
		for (unsigned short counter{0}; counter < 9; counter++) {
			feature_vector.push_back(0.0);
		}
		return;
	}

	// Add column maximum feature.
	feature_vector.push_back(maxima_(col));

	// Add column minimum feature.
	feature_vector.push_back(minima_(col));

	// Add column mean feature.
	feature_vector.push_back(means_(col));

	// Add column deviation feature.
	feature_vector.push_back(deviations_(col));

	// Add column uniqueness feature.
	double uniqueness{(matrix.col(col) - matrix(row, col)).abs().maxCoeff()};
	feature_vector.push_back(uniqueness);

	// Add column divergence feature.
	double divergence{(matrix.col(col) - matrix(row, col)).abs().sum()};
	if (((matrix.rows() - 1) * (matrix(row, col) - means_(col))) != 0) {
		divergence /= ((matrix.rows() - 1) * (matrix(row, col) - means_(col)));
	}
	feature_vector.push_back(divergence);

	// Add column leader feature.
	feature_vector.push_back(leaders_(col));

	// Add column interval feature.
	feature_vector.push_back(intervals_(col));

	// Add column outlierness feature.
	double outlierness{matrix(row, col)};
	if (means_(col) != deviations_(col)) {
		outlierness /= (means_(col) - deviations_(col));
	}
	feature_vector.push_back(outlierness);
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
populate_lsape_instance_(const GEDGraph & g, const GEDGraph & h, std::size_t num_threads) {
	DMatrix lsape_instance;
	lsape_method_->populate_instance(g, h, lsape_instance);
	lsape_instance_ = lsape_instance.matrix();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
compute_global_features_(const Eigen::ArrayXXd & substitution_matrix) {
	global_features_.clear();
	global_features_.push_back(substitution_matrix.maxCoeff());
	global_features_.push_back(substitution_matrix.minCoeff());
	global_features_.push_back(substitution_matrix.mean());
	if ((substitution_matrix.rows() * substitution_matrix.cols()) <= 1) {
		global_features_.push_back(0.0);
	}
	else {
		global_features_.push_back(std::sqrt((substitution_matrix - global_features_.at(2)).square().sum() / (substitution_matrix.rows() * substitution_matrix.cols() - 1)));
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
add_global_features_(std::vector<double> & feature_vector) const {
	for (auto feature = global_features_.begin(); feature != global_features_.end(); feature++) {
		feature_vector.push_back(*feature);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BipartiteML<UserNodeLabel, UserEdgeLabel>::
add_cell_features_(std::size_t row, std::size_t col, double node_cost, std::vector<double> & feature_vector) const {
	// Add node cost feature.
	feature_vector.push_back(node_cost);

	// Add edge cost feature.
	feature_vector.push_back(lsape_instance_(row, col) - node_cost);
}

}

#endif /* SRC_METHODS_BIPARTITE_ML_IPP_ */
