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
 * @file  ipfp.ipp
 * @brief ged::IPFP class definition.
 */

#ifndef SRC_METHODS_IPFP_IPP_
#define SRC_METHODS_IPFP_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
IPFP<UserNodeLabel, UserEdgeLabel>::
~IPFP() {}

template<class UserNodeLabel, class UserEdgeLabel>
IPFP<UserNodeLabel, UserEdgeLabel>::
IPFP(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
quadratic_model_{QAPE},
lsape_model_{LSAPESolver::ECBP},
epsilon_{0.001},
max_itrs_{100},
time_limit_in_sec_{0.0},
omega_{0.0},
qap_instance_() {}

// === Definitions of member functions inherited from LSBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) {

	// Start timer.
	Timer timer(time_limit_in_sec_);

	// Initialize output node map and upper bound.
	output_node_map = initial_node_map;
	double upper_bound{output_node_map.induced_cost()};

	// Initialize fractional solution x (is integral at start).
	DMatrix x;
	node_map_to_matrix_(initial_node_map, x);
	double linear_cost_x{compute_induced_linear_cost_(qap_instance_, x)};
	double overall_cost_x{initial_node_map.induced_cost()};
	bool x_is_integral{true};

	// Initialize gradient direction b.
	DMatrix b(qap_instance_.num_rows(), qap_instance_.num_cols());
	double linear_cost_b{0.0};
	double overall_cost_b{0.0};

	// Initialize alpha, beta, and the step width.
	double alpha{0.0};
	double beta{0.0};
	double step_width{0.0};

	// Initialize linear problem and solvers.
	DMatrix linear_problem(qap_instance_.num_rows(), qap_instance_.num_cols());
	LSAPESolver lsape_solver;
	lsape_solver.set_model(lsape_model_);
	LSAPSolver lsap_solver;
	if (quadratic_model_ == QAPE) {
		lsape_solver.set_problem(&linear_problem);
	}
	else {
		lsap_solver.set_problem(&linear_problem);
	}
	double min_linear_problem{0.0};

	// Initialize linear cost matrix.
	DMatrix linear_cost_matrix(qap_instance_.num_rows(), qap_instance_.num_cols());
	init_linear_cost_matrix_(qap_instance_, linear_cost_matrix);


	// Main loop.
	for (std::size_t current_itr{1}; not termination_criterion_met_(timer, alpha, min_linear_problem, current_itr, lower_bound, upper_bound); current_itr++) {


		// Compute the next gradient direction b and update the upper bound and the output node map.
		init_next_linear_problem_(qap_instance_, x, linear_cost_matrix, linear_problem);
		if (quadratic_model_ == QAPE) {
			solve_linear_problem_(qap_instance_, lsape_solver, min_linear_problem, linear_cost_b, overall_cost_b, b);
			if (overall_cost_b < upper_bound) {
				upper_bound = overall_cost_b;
				util::construct_node_map_from_solver(lsape_solver, output_node_map);
				this->ged_data_.compute_induced_cost(g, h, output_node_map);
				if (std::fabs(upper_bound - output_node_map.induced_cost()) > 0.000001) {
					throw Error("upper_bound: " + std::to_string(upper_bound) + ", linear_cost_b: " + std::to_string(linear_cost_b) + ", quadratic_cost_b: " + std::to_string(2*(upper_bound - linear_cost_b)) + ". Induced cost: " + std::to_string(output_node_map.induced_cost()) + ".");
				}
				output_node_map.set_induced_cost(upper_bound);
			}
		}
		else {
			solve_linear_problem_(qap_instance_, lsap_solver, min_linear_problem, linear_cost_b, overall_cost_b, b);
			if (overall_cost_b < upper_bound) {
				upper_bound = overall_cost_b;
				util::construct_node_map_from_solver(lsap_solver, output_node_map);
				output_node_map.set_induced_cost(upper_bound);
			}
		}

		// Determine the step width.
		alpha = min_linear_problem - 2 * overall_cost_x + linear_cost_x;
		beta = overall_cost_b + overall_cost_x - min_linear_problem - linear_cost_x;
		if (beta > 0.00001) {
			step_width = -alpha / (2 * beta);
		}

		// Update the possibly fractional solution x.
		if (beta <= 0.00001 or step_width >= 1) {
			x_is_integral = true;
			x = b;
			overall_cost_x = overall_cost_b;
			linear_cost_x = linear_cost_b;
		}
		else {
			x_is_integral = false;
			x += (b - x) * step_width;
			overall_cost_x -= (alpha * alpha) / (4 * beta);
			linear_cost_x = compute_induced_linear_cost_(qap_instance_, x);
		}
	}


	// If converged solution x is fractional, project it to integral solution.
	if (not x_is_integral) {
		DMatrix projection_problem(qap_instance_.num_rows(), qap_instance_.num_cols(), 1.0);
		projection_problem -= x;
		NodeMap projected_node_map(output_node_map);
		if (quadratic_model_ == QAPE) {
			lsape_solver.set_problem(&projection_problem);
			lsape_solver.solve();
			util::construct_node_map_from_solver(lsape_solver, projected_node_map);
		}
		else {
			lsap_solver.set_problem(&projection_problem);
			lsap_solver.solve();
			util::construct_node_map_from_solver(lsap_solver, projected_node_map);
		}
		this->ged_data_.compute_induced_cost(g, h, projected_node_map);
		if (projected_node_map.induced_cost() < output_node_map.induced_cost()) {
			output_node_map = projected_node_map;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
ls_runtime_init_(const GEDGraph & g, const GEDGraph & h) {
	omega_ = this->ged_data_.max_edit_cost(g, h) + 10.0;
	qap_instance_.init(this, g, h);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
ls_set_default_options_() {
	quadratic_model_ = QAPE;
	lsape_model_ = LSAPESolver::ECBP;
	epsilon_ = 0.001;
	max_itrs_ = 100;
	time_limit_in_sec_ = 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
IPFP<UserNodeLabel, UserEdgeLabel>::
ls_valid_options_string_() const {
	return "[--lsape-model <arg>] [--quadratic-model <arg>] [--iterations <arg>] [--time-limit <arg>] [--epsilon <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
IPFP<UserNodeLabel, UserEdgeLabel>::
ls_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "lsape-model") {
		if (arg == "EBP") {
			lsape_model_ = LSAPESolver::EBP;
		}
		else if (arg  == "FLWC") {
			lsape_model_ = LSAPESolver::FLWC;
		}
		else if (arg  == "FLCC") {
			lsape_model_ = LSAPESolver::FLCC;
		}
		else if (arg  == "FBP") {
			lsape_model_ = LSAPESolver::FBP;
		}
		else if (arg == "SFBP") {
			lsape_model_ = LSAPESolver::SFBP;
		}
		else if (arg  == "FBP0") {
			lsape_model_ = LSAPESolver::FBP0;
		}
		else if (arg  == "ECBP") {
			lsape_model_ = LSAPESolver::ECBP;
		}
		else {
			throw ged::Error(std::string("Invalid argument ") + arg + " for option lsape-model. Usage: options = \"[--lsape-model ECBP|EBP|FLWC|FLCC|FBP|SFBP|FBP0] [...]\"");
		}
		return true;
	}
	else if (option == "time-limit") {
		try {
			time_limit_in_sec_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option time-limit.  Usage: options = \"[--time-limit <convertible to double>] [...]");
		}
		return true;
	}
	else if (option == "quadratic-model") {
		if (arg == "QAPE") {
			quadratic_model_ = QAPE;
		}
		else if (arg  == "B-QAP") {
			quadratic_model_ = B_QAP;
		}
		else if (arg  == "C-QAP") {
			quadratic_model_ = C_QAP;
		}
		else {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option quadratic-model. Usage: options = \"[--quadratic-model QAPE|B-QAP|C-QAP] [...]\"");
		}
		return true;
	}
	else if (option == "iterations") {
		try {
			max_itrs_ = std::stoi(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option iterations. Usage: options = \"[--iterations <convertible to int>] [...]");
		}
		return true;
	}
	else if (option == "epsilon") {
		try {
			epsilon_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option epsilon. Usage: options = \"[--epsilon <convertible to double greater 0>] [...]");
		}
		if (epsilon_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option epsilon. Usage: options = \"[--epsilon <convertible to double greater 0>] [...]");
		}
		return true;
	}
	return false;
}

// == Definitions of private helper member functions. ===

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
node_map_to_matrix_(const NodeMap & node_map, DMatrix & matrix) const {
	matrix.resize(qap_instance_.num_rows(), qap_instance_.num_cols());
	matrix.set_to_val(0.0);
	std::vector<NodeMap::Assignment> relation;
	node_map.as_relation(relation);
	for (auto assignment : relation) {
		std::size_t row{undefined()};
		std::size_t col{undefined()};
		if (assignment.first != GEDGraph::dummy_node()) {
			row = assignment.first;
		}
		if (assignment.second != GEDGraph::dummy_node()) {
			col = assignment.second;
		}
		if ((row != undefined()) and (col != undefined())) {
			matrix(row, col) = 1.0;
		}
		else if (row != undefined()) {
			if (quadratic_model_ == QAPE) {
				matrix(row, matrix.num_cols() - 1) = 1.0;
			}
			else if (quadratic_model_ == B_QAP){
				matrix(row, row + qap_instance_.num_nodes_h()) = 1.0;
			}
		}
		else if (col != undefined()) {
			if (quadratic_model_ == QAPE) {
				matrix(matrix.num_rows() - 1, col) = 1.0;
			}
			else if (quadratic_model_ == B_QAP){
				matrix(col + qap_instance_.num_nodes_g(), col) = 1.0;
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
init_linear_cost_matrix_(const QAPInstance_ & qap_instance, DMatrix & linear_cost_matrix) const {
	for (std::size_t row_1{0}; row_1 < linear_cost_matrix.num_rows(); row_1++) {
		for (std::size_t col_1{0}; col_1 < linear_cost_matrix.num_cols(); col_1++) {
			linear_cost_matrix(row_1, col_1) = qap_instance(row_1, col_1);
		}
	}
}


template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
init_next_linear_problem_(const QAPInstance_ & qap_instance, const DMatrix & x, const DMatrix & linear_cost_matrix, DMatrix & linear_problem) const {
	std::vector<NodeMap::Assignment> support_x;
	for (std::size_t col{0}; col < linear_problem.num_cols(); col++) {
		for (std::size_t row{0}; row < linear_problem.num_rows(); row++) {
			if (x(row, col) != 0) {
				support_x.emplace_back(row, col);
			}
		}
	}
	for (std::size_t col{0}; col < linear_problem.num_cols(); col++) {
		for (std::size_t row{0}; row < linear_problem.num_rows(); row++) {
			linear_problem(row, col) = 0.0;
			for (const auto & assignment : support_x) {
				linear_problem(row, col) += qap_instance(row, col, assignment.first, assignment.second) * x(assignment.first, assignment.second);
			}
		}
	}
	linear_problem += linear_cost_matrix;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
IPFP<UserNodeLabel, UserEdgeLabel>::
termination_criterion_met_(const Timer & timer, const double & alpha, const double & min_linear_problem, const std::size_t & current_itr, double lower_bound, double upper_bound) const {
	if (current_itr == 1) {
		return false;
	}
	if (lower_bound >= upper_bound) {
		return true;
	}
	if (timer.expired()) {
		return true;
	}
	if (max_itrs_ > 0 and current_itr > max_itrs_) {
		return true;
	}
	if (min_linear_problem < 0.00001) {
		return std::fabs(alpha) < epsilon_;
	}
	return std::fabs(alpha / min_linear_problem) < epsilon_;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
solve_linear_problem_(const QAPInstance_ & qap_instance, LSAPSolver & solver,
		double & min_linear_problem, double & linear_cost_b, double & overall_cost_b, DMatrix & b) const {
	solver.solve();
	min_linear_problem = solver.minimal_cost();
	solver_to_matrix_(solver, b);
	linear_cost_b = compute_induced_linear_cost_(qap_instance, solver);
	overall_cost_b = linear_cost_b + 0.5 * compute_induced_quadratic_cost_(qap_instance, solver);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
solve_linear_problem_(const QAPInstance_ & qap_instance, LSAPESolver & solver,
		double & min_linear_problem, double & linear_cost_b, double & overall_cost_b, DMatrix & b) const {
	solver.solve();
	min_linear_problem = solver.minimal_cost();
	solver_to_matrix_(solver, b);
	linear_cost_b = compute_induced_linear_cost_(qap_instance, solver);
	overall_cost_b = linear_cost_b + 0.5 * compute_induced_quadratic_cost_(qap_instance, solver);
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
compute_induced_linear_cost_(const QAPInstance_ & qap_instance, const DMatrix & x) const {
	double linear_cost{0.0};
	for (std::size_t row{0}; row < qap_instance.num_rows(); row++) {
		for (std::size_t col{0}; col < qap_instance.num_cols(); col++) {
			linear_cost += qap_instance(row, col) * x(row, col);
		}
	}
	return linear_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
solver_to_matrix_(const LSAPSolver & solver, DMatrix & b) const {
	b.matrix().setZero();
	for (std::size_t row{0}; row < solver.num_rows(); row++) {
		if (solver.get_assigned_col(row) < solver.num_cols()) {
			b(row, solver.get_assigned_col(row)) = 1.0;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
solver_to_matrix_(const LSAPESolver & solver, DMatrix & b) const {
	b.matrix().setZero();
	for (std::size_t row{0}; row < solver.num_rows(); row++) {
		b(row, solver.get_assigned_col(row)) = 1.0;
	}
	for (std::size_t col{0}; col < solver.num_cols(); col++) {
		if (solver.get_assigned_row(col) == solver.num_rows()) {
			b(solver.get_assigned_row(col), col) = 1.0;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
compute_induced_linear_cost_(const QAPInstance_ & qap_instance, const LSAPSolver & solver) const {
	double linear_cost{0.0};
	for (std::size_t row{0}; row < solver.num_rows(); row++) {
		if (solver.get_assigned_col(row) < solver.num_cols()) {
			linear_cost += qap_instance(row, solver.get_assigned_col(row));
		}
	}
	return linear_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
compute_induced_linear_cost_(const QAPInstance_ & qap_instance, const LSAPESolver & solver) const {
	double linear_cost{0.0};
	for (std::size_t row{0}; row < solver.num_rows(); row++) {
		linear_cost += qap_instance(row, solver.get_assigned_col(row));
	}
	for (std::size_t col{0}; col < solver.num_cols(); col++) {
		if (solver.get_assigned_row(col) == solver.num_rows()) {
			linear_cost += qap_instance(solver.get_assigned_row(col), col);
		}
	}
	return linear_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
compute_induced_quadratic_cost_(const QAPInstance_ & qap_instance, const LSAPSolver & solver) const {
	double quadratic_cost{0.0};
	std::vector<NodeMap::Assignment> assignments;
	for (std::size_t row{0}; row < solver.num_rows(); row++) {
		if (solver.get_assigned_col(row) < solver.num_cols()) {
			assignments.emplace_back(row, solver.get_assigned_col(row));
		}
	}
	for (const auto & assignment_1 : assignments) {
		for (const auto & assignment_2 : assignments) {
			quadratic_cost += qap_instance(assignment_1.first, assignment_1.second, assignment_2.first, assignment_2.second);
		}
	}
	return quadratic_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
compute_induced_quadratic_cost_(const QAPInstance_ & qap_instance, const LSAPESolver & solver) const {
	double quadratic_cost{0.0};
	std::vector<NodeMap::Assignment> assignments;
	for (std::size_t row{0}; row < solver.num_rows(); row++) {
		assignments.emplace_back(row, solver.get_assigned_col(row));
	}
	for (std::size_t col{0}; col < solver.num_cols(); col++) {
		if (solver.get_assigned_row(col) == solver.num_rows()) {
			assignments.emplace_back(solver.get_assigned_row(col), col);
		}
	}
	for (const auto & assignment_1 : assignments) {
		for (const auto & assignment_2 : assignments) {
			quadratic_cost += qap_instance(assignment_1.first, assignment_1.second, assignment_2.first, assignment_2.second);
		}
	}
	return quadratic_cost;
}

// === Definition of private class QAPInstance_. ===
template<class UserNodeLabel, class UserEdgeLabel>
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
QAPInstance_() :
ipfp_{nullptr},
g_{nullptr},
h_{nullptr},
num_nodes_g_{0},
num_nodes_h_{0},
translation_factor_{0.0} {}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
init(const IPFP<UserNodeLabel, UserEdgeLabel> * ipfp, const GEDGraph & g, const GEDGraph & h)  {
	ipfp_ = ipfp;
	g_ = &g;
	h_ = &h;
	num_nodes_g_ = g_->num_nodes();
	num_nodes_h_ = h_->num_nodes();
	if (ipfp_->quadratic_model_ == C_QAP) {
		if (num_nodes_g_ > num_nodes_h_) {
			translation_factor_ = 3 * std::max(ipfp_->ged_data_.max_node_del_cost(*g_), ipfp_->ged_data_.max_edge_del_cost(*g_));
		}
		else if (num_nodes_g_ < num_nodes_h_) {
			translation_factor_ = 3 * std::max(ipfp_->ged_data_.max_node_ins_cost(*h_), ipfp_->ged_data_.max_edge_ins_cost(*h_));
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
num_rows() const {
	switch(ipfp_->quadratic_model_) {
	case QAPE:
		return static_cast<std::size_t>(num_nodes_g_) + 1;
	case B_QAP:
		return static_cast<std::size_t>(num_nodes_g_) + static_cast<std::size_t>(num_nodes_h_);
	default:
		return static_cast<std::size_t>(num_nodes_g_);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
num_cols() const {
	switch(ipfp_->quadratic_model_) {
	case QAPE:
		return static_cast<std::size_t>(num_nodes_h_) + 1;
	case B_QAP:
		return static_cast<std::size_t>(num_nodes_h_) + static_cast<std::size_t>(num_nodes_g_);
	default:
		return static_cast<std::size_t>(num_nodes_h_);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
num_nodes_g() const {
	return num_nodes_g_;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
num_nodes_h() const {
	return num_nodes_h_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
operator() (std::size_t row, std::size_t col) const {
	if (row >= num_rows() or col >= num_cols()) {
		throw Error("Out of range error for FrankWolfe<UserNodeLabel, UserEdgeLabel>::QAPInstance_.");
	}
	double return_val{0.0};
	switch(ipfp_->quadratic_model_) {
	case QAP:
		return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), h_->get_node_label(col));
		break;
	case C_QAP:
		return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), h_->get_node_label(col));
		if (num_nodes_g_ > num_nodes_h_) {
			return_val -= ipfp_->ged_data_.node_cost(g_->get_node_label(row), dummy_label());
		}
		else if (num_nodes_g_ < num_nodes_h_) {
			return_val -= ipfp_->ged_data_.node_cost(dummy_label(), h_->get_node_label(col));
		}
		if (num_nodes_g_ != num_nodes_h_) {
			return_val += translation_factor_;
		}
		break;
	case QAPE:
		if (row < num_nodes_g_ and col < num_nodes_h_) {
			return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), h_->get_node_label(col));
		}
		else if (row < num_nodes_g_) {
			return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), dummy_label());
		}
		else if (col < num_nodes_h_) {
			return_val += ipfp_->ged_data_.node_cost(dummy_label(), h_->get_node_label(col));
		}
		break;
	case B_QAP:
		if (row < num_nodes_g_ and col >= num_nodes_h_ and col != row + num_nodes_h_) {
			return_val += ipfp_->omega_;
		}
		else if (row >= num_nodes_g_ and col < num_nodes_h_ and row != col + num_nodes_g_) {
			return_val += ipfp_->omega_;
		}
		else if (row < num_nodes_g_ and col < num_nodes_h_) {
			return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), h_->get_node_label(col));
		}
		else if (row < num_nodes_g_) {
			return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), dummy_label());
		}
		else if (col < num_nodes_h_) {
			return_val += ipfp_->ged_data_.node_cost(dummy_label(), h_->get_node_label(col));
		}
		break;
	}
	return return_val;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
quadratic_cost_b_qap_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const {
	if (row_1 == row_2 or col_1 == col_2) {
		return 0.0;
	}
	double return_val{0.0};
	if (row_1 < num_nodes_g_ and  col_1 >= num_nodes_h_ and col_1 != row_1 + num_nodes_h_) {
		return_val += ipfp_->omega_;
	}
	else if (row_2 < num_nodes_g_ and col_2 >= num_nodes_h_ and col_2 != row_2 + num_nodes_h_) {
		return_val += ipfp_->omega_;
	}
	else if (col_1 < num_nodes_h_ and row_1 >= num_nodes_g_ and row_1 != col_1 + num_nodes_g_) {
		return_val += ipfp_->omega_;
	}
	else if (col_2 < num_nodes_h_ and row_2 >= num_nodes_g_ and row_2 != col_2 + num_nodes_g_) {
		return_val += ipfp_->omega_;
	}
	else if (row_1 < num_nodes_g_ and col_1 < num_nodes_h_ and row_2 < num_nodes_g_ and col_2 < num_nodes_h_) {
		GEDGraph::EdgeID edge_g(g_->get_edge(row_1, row_2));
		GEDGraph::EdgeID edge_h(h_->get_edge(col_1, col_2));
		if (edge_g != GEDGraph::dummy_edge() and edge_h != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), h_->get_edge_label(edge_h));
		}
		else if (edge_g != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), dummy_label());
		}
		else if (edge_h != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(edge_h));
		}
	}
	else if (row_1 < num_nodes_g_ and row_2 < num_nodes_g_) {
		GEDGraph::EdgeID edge_g(g_->get_edge(row_1, row_2));
		if (edge_g != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), dummy_label());
		}
	}
	else if (col_1 < num_nodes_h_ and col_2 < num_nodes_h_) {
		GEDGraph::EdgeID edge_h(h_->get_edge(col_1, col_2));
		if (edge_h != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(edge_h));
		}
	}
	return return_val;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
quadratic_cost_c_qap_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const {
	if (row_1 == row_2 or col_1 == col_2) {
		return 0.0;
	}
	double return_val{0.0};
	GEDGraph::EdgeID edge_g(g_->get_edge(row_1, row_2));
	GEDGraph::EdgeID edge_h(h_->get_edge(col_1, col_2));
	if (edge_g != GEDGraph::dummy_edge() and edge_h != GEDGraph::dummy_edge()) {
		return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), h_->get_edge_label(edge_h));
	}
	else if (edge_g != GEDGraph::dummy_edge()) {
		return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), dummy_label());
	}
	else if (edge_h != GEDGraph::dummy_edge()) {
		return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(edge_h));
	}
	if (edge_g != GEDGraph::dummy_edge() and num_nodes_g_ > num_nodes_h_) {
		return_val -= 3 * ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), dummy_label());
	}
	else if (edge_h != GEDGraph::dummy_edge() and num_nodes_g_ < num_nodes_h_) {
		return_val -= 3 * ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(edge_h));
	}
	if (num_nodes_g_ != num_nodes_h_) {
		return_val += translation_factor_;
	}
	return return_val;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
quadratic_cost_qap_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const {
	if (row_1 == row_2 or col_1 == col_2) {
		return 0.0;
	}
	double return_val{0.0};
	GEDGraph::EdgeID edge_g(g_->get_edge(row_1, row_2));
	GEDGraph::EdgeID edge_h(h_->get_edge(col_1, col_2));
	if (edge_g != GEDGraph::dummy_edge() and edge_h != GEDGraph::dummy_edge()) {
		return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), h_->get_edge_label(edge_h));
	}
	else if (edge_g != GEDGraph::dummy_edge()) {
		return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), dummy_label());
	}
	else if (edge_h != GEDGraph::dummy_edge()) {
		return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(edge_h));
	}
	return return_val;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
quadratic_cost_qape_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const {
	//std::cout << "1\n"; // todo rm
	if (row_1 == row_2 and col_1 == col_2) {
		return 0.0;
	}
	if (row_1 == row_2 and row_1 < num_nodes_g_) {
		return 0.0;
	}
	if (col_1 == col_2 and col_1 < num_nodes_h_) {
		return 0.0;
	}
	double return_val{0.0};
	if (row_1 < num_nodes_g_ and col_1 < num_nodes_h_ and row_2 < num_nodes_g_ and col_2 < num_nodes_h_) {
		GEDGraph::EdgeID edge_g(g_->get_edge(row_1, row_2));
		GEDGraph::EdgeID edge_h(h_->get_edge(col_1, col_2));
		if (edge_g != GEDGraph::dummy_edge() and edge_h != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), h_->get_edge_label(edge_h));
		}
		else if (edge_g != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), dummy_label());
		}
		else if (edge_h != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(edge_h));
		}
	}
	else if (row_1 < num_nodes_g_ and row_2 < num_nodes_g_) {
		GEDGraph::EdgeID edge_g(g_->get_edge(row_1, row_2));
		if (edge_g != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(edge_g), dummy_label());
		}
	}
	else if (col_1 < num_nodes_h_ and col_2 < num_nodes_h_) {
		GEDGraph::EdgeID edge_h(h_->get_edge(col_1, col_2));
		if (edge_h != GEDGraph::dummy_edge()) {
			return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(edge_h));
		}
	}
	//std::cout << "2\n"; // todo rm
	return return_val;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
operator() (std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const {
	if (row_1 >= num_rows() or col_1 >= num_cols() or row_2 >= num_rows() or col_2 >= num_cols()) {
		throw Error("Out of range error for FrankWolfe<UserNodeLabel, UserEdgeLabel>::QAPInstance_.");
	}
	double return_val{0.0};
	switch(ipfp_->quadratic_model_) {
	case QAP:
		return_val += quadratic_cost_qap_(row_1, col_1, row_2, col_2);
		break;
	case C_QAP:
		return_val += quadratic_cost_c_qap_(row_1, col_1, row_2, col_2);
		break;
	case QAPE:
		return_val += quadratic_cost_qape_(row_1, col_1, row_2, col_2);
		break;
	case B_QAP:
		return_val += quadratic_cost_b_qap_(row_1, col_1, row_2, col_2);
		break;
	}
	return return_val;
}

}

#endif /* SRC_METHODS_IPFP_IPP_ */
