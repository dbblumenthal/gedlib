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
epsilon_{0.000001},
max_itrs_{100},
time_limit_in_sec_{0.0},
omega_{0.0},
qap_instance_() {}

// === Definitions of member functions inherited from LSBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) {

	// Initialize variables.
	DMatrix x_0;
	node_map_to_matrix_(initial_node_map, x_0);
	Timer timer(time_limit_in_sec_);
	DMatrix x(x_0);
	DMatrix b(x_0);
	double linear_cost_x{compute_induced_linear_cost_(qap_instance_, x)};
	double linear_cost_b{0.0};
	double overall_cost_x{linear_cost_x + 0.5 * compute_induced_quadratic_cost_(qap_instance_, x)};
	double overall_cost_b{0.0};
	double min_linear_problem{0.0};
	bool x_is_integral{true};
	double alpha{0.0};
	double beta{0.0};
	double step_width{0.0};
	double upper_bound{initial_node_map.induced_cost()};
	output_node_map = initial_node_map;

	// Initialize linear problem and solvers.
	DMatrix linear_problem(qap_instance_.num_rows(), qap_instance_.num_cols());
	LSAPESolver lsape_solver;
	lsape_solver.set_model(lsape_model_);
	LSAPSolver lsap_solver;

	// Main loop.
	for (std::size_t current_itr{1}; not termination_criterion_met_(timer, alpha, min_linear_problem, current_itr, lower_bound, upper_bound); current_itr++) {
		init_next_linear_problem_(qap_instance_, x, linear_problem);
		if (quadratic_model_ == QAPE) {
			solve_linear_problem_(linear_problem, qap_instance_, lsape_solver, min_linear_problem, linear_cost_b, overall_cost_b);
			if (overall_cost_b < upper_bound) {
				upper_bound = overall_cost_b;
				util::construct_node_map_from_solver(lsape_solver, output_node_map);
			}
		}
		else {
			solve_linear_problem_(linear_problem, qap_instance_, lsap_solver, min_linear_problem, linear_cost_b, overall_cost_b);
			if (overall_cost_b < upper_bound) {
				upper_bound = overall_cost_b;
				util::construct_node_map_from_solver(lsap_solver, output_node_map);
			}
		}
		alpha = min_linear_problem - 2 * overall_cost_x + linear_cost_x;
		beta = overall_cost_b + overall_cost_x - min_linear_problem - linear_cost_x;
		if (beta > 0.00001) {
			step_width = -alpha / (2 * beta);
		}
		if (beta <= 0.00001 or step_width >= 1) {
			x_is_integral = true;
			if (quadratic_model_ == QAPE) {
				solver_to_matrix_(lsape_solver, x);
			}
			else {
				solver_to_matrix_(lsap_solver, x);
			}
			overall_cost_x = overall_cost_b;
			linear_cost_x = linear_cost_b;
		}
		else {
			x_is_integral = false;
			if (quadratic_model_ == QAPE) {
				solver_to_matrix_(lsape_solver, b);
			}
			else {
				solver_to_matrix_(lsap_solver, b);
			}
			x += (b - x) * step_width;
			overall_cost_x -= (alpha * alpha) / (4 * beta);
			linear_cost_x = compute_induced_linear_cost_(qap_instance_, x);
		}
	}

	// Project doubly stochastic matrix to integral solution.
	if (not x_is_integral) {
		DMatrix projection_problem(qap_instance_.num_rows(), qap_instance_.num_cols(), 1.0);
		projection_problem -= x;
		if (quadratic_model_ == QAPE) {
			solve_linear_problem_(projection_problem, qap_instance_, lsape_solver, min_linear_problem, linear_cost_b, overall_cost_b);
			if (overall_cost_b < upper_bound) {
				upper_bound = overall_cost_b;
				util::construct_node_map_from_solver(lsape_solver, output_node_map);
			}
		}
		else {
			solve_linear_problem_(projection_problem, qap_instance_, lsap_solver, min_linear_problem, linear_cost_b, overall_cost_b);
			if (overall_cost_b < upper_bound) {
				upper_bound = overall_cost_b;
				util::construct_node_map_from_solver(lsap_solver, output_node_map);
			}
		}
	}

	// Set the induced cost of the output node map to the computed upper bound.
	output_node_map.set_induced_cost(upper_bound);
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
init_next_linear_problem_(const QAPInstance_ & qap_instance, const DMatrix & x_k, DMatrix & linear_problem) const {
	for (std::size_t row_1{0}; row_1 < linear_problem.num_rows(); row_1++) {
		for (std::size_t col_1{0}; col_1 < linear_problem.num_cols(); col_1++) {
			linear_problem(row_1, col_1) = 0.0;
			if (x_k(row_1, col_1) > 0.0) {
				for (std::size_t row_2{0}; row_2 < linear_problem.num_rows(); row_2++) {
					for (std::size_t col_2{0}; col_2 < linear_problem.num_cols(); col_2++) {
						linear_problem(row_1, col_1) += qap_instance(row_1, col_1, row_2, col_2);
					}
				}
				linear_problem(row_1, col_1) *= x_k(row_1, col_1);
			}
		}
	}
	for (std::size_t row_1{0}; row_1 < linear_problem.num_rows(); row_1++) {
		for (std::size_t col_1{0}; col_1 < linear_problem.num_cols(); col_1++) {
			linear_problem(row_1, col_1) += qap_instance(row_1, col_1);
		}
	}
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
solve_linear_problem_(const DMatrix & linear_problem, const QAPInstance_ & qap_instance, LSAPSolver & solver,
		double & min_linear_problem, double & linear_cost_b, double & overall_cost_b) const {
	solver.set_problem(linear_problem);
	solver.solve();
	min_linear_problem = solver.minimal_cost();
	linear_cost_b = compute_induced_linear_cost_(qap_instance, solver);
	overall_cost_b = linear_cost_b + 0.5 * compute_induced_quadratic_cost_(qap_instance, solver);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
solve_linear_problem_(const DMatrix & linear_problem, const QAPInstance_ & qap_instance, LSAPESolver & solver,
		double & min_linear_problem, double & linear_cost_b, double & overall_cost_b) const {
	solver.set_problem(linear_problem);
	solver.solve();
	min_linear_problem = solver.minimal_cost();
	linear_cost_b = compute_induced_linear_cost_(qap_instance, solver);
	overall_cost_b = linear_cost_b + 0.5 * compute_induced_quadratic_cost_(qap_instance, solver);
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
compute_induced_linear_cost_(const QAPInstance_ & qap_instance, const DMatrix & x_k) const {
	double linear_cost{0.0};
	for (std::size_t row{0}; row < qap_instance.num_rows(); row++) {
		for (std::size_t col{0}; col < qap_instance.num_cols(); col++) {
			linear_cost += qap_instance(row, col) * x_k(row, col);
		}
	}
	return linear_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
solver_to_matrix_(const LSAPSolver & solver, DMatrix & b_k_plus_1) const {
	for (std::size_t row{0}; row < b_k_plus_1.num_rows(); row++) {
		for (std::size_t col{0}; col < b_k_plus_1.num_cols(); col++) {
			b_k_plus_1(row, col) = 0.0;
		}
	}
	for (std::size_t row{0}; row < solver.num_rows(); row++) {
		if (solver.get_assigned_col(row) < solver.num_cols()) {
			b_k_plus_1(row, solver.get_assigned_col(row)) = 1.0;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
solver_to_matrix_(const LSAPESolver & solver, DMatrix & b_k_plus_1) const {
	for (std::size_t row{0}; row < b_k_plus_1.num_rows(); row++) {
		for (std::size_t col{0}; col < b_k_plus_1.num_cols(); col++) {
			b_k_plus_1(row, col) = 0.0;
		}
	}
	for (std::size_t row{0}; row < solver.num_rows(); row++) {
		b_k_plus_1(row, solver.get_assigned_col(row)) = 1.0;
	}
	for (std::size_t col{0}; col < solver.num_cols(); col++) {
		if (solver.get_assigned_row(col) == solver.num_rows()) {
			b_k_plus_1(solver.get_assigned_row(col), col) = 1.0;
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
compute_induced_quadratic_cost_(const QAPInstance_ & qap_instance, const DMatrix & x_k) const {
	double quadratic_cost{0.0};
	for (std::size_t row_1{0}; row_1 < qap_instance.num_rows(); row_1++) {
		for (std::size_t col_1{0}; col_1 < qap_instance.num_cols(); col_1++) {
			if (x_k(row_1, col_1) == 0.0) {
				continue;
			}
			for (std::size_t row_2{0}; row_2 < qap_instance.num_rows(); row_2++) {
				for (std::size_t col_2{0}; col_2 < qap_instance.num_cols(); col_2++) {
					quadratic_cost += qap_instance(row_1, col_1, row_2, col_2) * x_k(row_2, col_2);
				}
			}
		}
	}
	return quadratic_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
compute_induced_quadratic_cost_(const QAPInstance_ & qap_instance, const LSAPSolver & solver) const {
	double quadratic_cost{0.0};
	for (std::size_t row_1{0}; row_1 < solver.num_rows(); row_1++) {
		if (solver.get_assigned_col(row_1) == solver.num_cols()) {
			continue;
		}
		for (std::size_t row_2{0}; row_2 < solver.num_rows(); row_2++) {
			if (solver.get_assigned_col(row_2) < solver.num_cols()) {
				quadratic_cost += qap_instance(row_1, solver.get_assigned_col(row_1), row_2, solver.get_assigned_col(row_2));
			}
		}
	}
	return quadratic_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
compute_induced_quadratic_cost_(const QAPInstance_ & qap_instance, const LSAPESolver & solver) const {
	double quadratic_cost{0.0};
	for (std::size_t row_1{0}; row_1 < solver.num_rows(); row_1++) {
		for (std::size_t row_2{0}; row_2 < solver.num_rows(); row_2++) {
			quadratic_cost += qap_instance(row_1, solver.get_assigned_col(row_1), row_2, solver.get_assigned_col(row_2));
		}
	}
	for (std::size_t col_1{0}; col_1 < solver.num_cols(); col_1++) {
		if (solver.get_assigned_row(col_1) == solver.num_rows()) {
			for (std::size_t row_2{0}; row_2 < solver.num_rows(); row_2++) {
				quadratic_cost += 2 * qap_instance(solver.get_assigned_row(col_1), col_1, row_2, solver.get_assigned_col(row_2));
			}
			for (std::size_t col_2{0}; col_2 < solver.num_cols(); col_2++) {
				if (solver.get_assigned_row(col_2) == solver.num_rows()) {
					quadratic_cost += 2 * qap_instance(solver.get_assigned_row(col_1), col_1, solver.get_assigned_row(col_2), col_2);
				}
			}
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
translation_factor_{0.0} {}

template<class UserNodeLabel, class UserEdgeLabel>
void
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
init(const IPFP<UserNodeLabel, UserEdgeLabel> * ipfp, const GEDGraph & g, const GEDGraph & h)  {
	ipfp_ = ipfp;
	g_ = &g;
	h_ = &h;
	if (ipfp_->quadratic_model_ == C_QAP) {
		if (g_->num_nodes() > h_->num_nodes()) {
			translation_factor_ = 3 * std::max(ipfp_->ged_data_.max_node_del_cost(*g_), ipfp_->ged_data_.max_edge_del_cost(*g_));
		}
		else if (g_->num_nodes() < h_->num_nodes()) {
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
		return static_cast<std::size_t>(g_->num_nodes()) + 1;
	case B_QAP:
		return static_cast<std::size_t>(g_->num_nodes()) + static_cast<std::size_t>(h_->num_nodes());
	default:
		return static_cast<std::size_t>(g_->num_nodes());
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
num_cols() const {
	switch(ipfp_->quadratic_model_) {
	case QAPE:
		return static_cast<std::size_t>(h_->num_nodes()) + 1;
	case B_QAP:
		return static_cast<std::size_t>(h_->num_nodes()) + static_cast<std::size_t>(g_->num_nodes());
	default:
		return static_cast<std::size_t>(h_->num_nodes());
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
num_nodes_g() const {
	return g_->num_nodes();
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
num_nodes_h() const {
	return h_->num_nodes();
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
		if (g_->num_nodes() > h_->num_nodes()) {
			return_val -= ipfp_->ged_data_.node_cost(g_->get_node_label(row), dummy_label());
		}
		else if (g_->num_nodes() < h_->num_nodes()) {
			return_val -= ipfp_->ged_data_.node_cost(dummy_label(), h_->get_node_label(col));
		}
		if (g_->num_nodes() != h_->num_nodes()) {
			return_val += translation_factor_;
		}
		break;
	case QAPE:
		if (row < g_->num_nodes() and col < h_->num_nodes()) {
			return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), h_->get_node_label(col));
		}
		else if (row < g_->num_nodes()) {
			return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), dummy_label());
		}
		else if (col < h_->num_nodes()) {
			return_val += ipfp_->ged_data_.node_cost(dummy_label(), h_->get_node_label(col));
		}
		break;
	case B_QAP:
		if (row < g_->num_nodes() and col >= h_->num_nodes() and col != row + h_->num_nodes()) {
			return_val += ipfp_->omega_;
		}
		else if (row >= g_->num_nodes() and col < h_->num_nodes() and row != col + g_->num_nodes()) {
			return_val += ipfp_->omega_;
		}
		else if (row < g_->num_nodes() and col < h_->num_nodes()) {
			return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), h_->get_node_label(col));
		}
		else if (row < g_->num_nodes()) {
			return_val += ipfp_->ged_data_.node_cost(g_->get_node_label(row), dummy_label());
		}
		else if (col < h_->num_nodes()) {
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
	double return_val{0.0};
	if (row_1 < g_->num_nodes() and  col_1 >= h_->num_nodes() and col_1 != row_1 + h_->num_nodes()) {
		return_val += ipfp_->omega_;
	}
	else if (row_2 < g_->num_nodes() and col_2 >= h_->num_nodes() and col_2 != row_2 + h_->num_nodes()) {
		return_val += ipfp_->omega_;
	}
	else if (col_1 < h_->num_nodes() and row_1 >= g_->num_nodes() and row_1 != col_1 + g_->num_nodes()) {
		return_val += ipfp_->omega_;
	}
	else if (col_2 < h_->num_nodes() and row_2 >= g_->num_nodes() and row_2 != col_2 + g_->num_nodes()) {
		return_val += ipfp_->omega_;
	}
	else if (row_1 < g_->num_nodes() and col_1 < h_->num_nodes() and row_2 < g_->num_nodes() and col_2 < h_->num_nodes()) {
		if (g_->is_edge(row_1, row_2) and h_->is_edge(col_1, col_2)) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), h_->get_edge_label(h_->get_edge(col_1, col_2)));
		}
		else if (g_->is_edge(row_1, row_2)) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), dummy_label());
		}
		else if (h_->is_edge(col_1, col_2)) {
			return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(h_->get_edge(col_1, col_2)));
		}
	}
	else if (row_1 < g_->num_nodes() and row_2 < g_->num_nodes()) {
		if (g_->is_edge(row_1, row_2)) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), dummy_label());
		}
	}
	else if (col_1 < h_->num_nodes() and col_2 < h_->num_nodes()) {
		if (h_->is_edge(col_1, col_2)) {
			return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(h_->get_edge(col_1, col_2)));
		}
	}
	return return_val;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
quadratic_cost_c_qap_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const {
	double return_val{0.0};
	if (g_->is_edge(row_1, row_2) and h_->is_edge(col_1, col_2)) {
		return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), h_->get_edge_label(h_->get_edge(col_1, col_2)));
	}
	else if (g_->is_edge(row_1, row_2)) {
		return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), dummy_label());
	}
	else if (h_->is_edge(col_1, col_2)) {
		return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(h_->get_edge(col_1, col_2)));
	}
	if (g_->is_edge(row_1, row_2) and g_->num_nodes() > h_->num_nodes()) {
		return_val -= 3 * ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), dummy_label());
	}
	else if (h_->is_edge(col_1, col_2) and g_->num_nodes() < h_->num_nodes()) {
		return_val -= 3 * ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(h_->get_edge(col_1, col_2)));
	}
	if (g_->num_nodes() != h_->num_nodes()) {
		return_val += translation_factor_;
	}
	return return_val;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
quadratic_cost_qap_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const {
	double return_val{0.0};
	if (g_->is_edge(row_1, row_2) and h_->is_edge(col_1, col_2)) {
		return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), h_->get_edge_label(h_->get_edge(col_1, col_2)));
	}
	else if (g_->is_edge(row_1, row_2)) {
		return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), dummy_label());
	}
	else if (h_->is_edge(col_1, col_2)) {
		return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(h_->get_edge(col_1, col_2)));
	}
	return return_val;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
IPFP<UserNodeLabel, UserEdgeLabel>::
QAPInstance_ ::
quadratic_cost_qape_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const {
	double return_val{0.0};
	if (row_1 < g_->num_nodes() and col_1 < h_->num_nodes() and row_2 < g_->num_nodes() and col_2 < h_->num_nodes()) {
		if (g_->is_edge(row_1, row_2) and h_->is_edge(col_1, col_2)) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), h_->get_edge_label(h_->get_edge(col_1, col_2)));
		}
		else if (g_->is_edge(row_1, row_2)) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), dummy_label());
		}
		else if (h_->is_edge(col_1, col_2)) {
			return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(h_->get_edge(col_1, col_2)));
		}
	}
	else if (row_1 < g_->num_nodes() and row_2 < g_->num_nodes()) {
		if (g_->is_edge(row_1, row_2)) {
			return_val += ipfp_->ged_data_.edge_cost(g_->get_edge_label(g_->get_edge(row_1, row_2)), dummy_label());
		}
	}
	else if (col_1 < h_->num_nodes() and col_2 < h_->num_nodes()) {
		if (h_->is_edge(col_1, col_2)) {
			return_val += ipfp_->ged_data_.edge_cost(dummy_label(), h_->get_edge_label(h_->get_edge(col_1, col_2)));
		}
	}
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
	if (row_1 == row_2 or col_1 == col_2) {
		return 0.0;
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
