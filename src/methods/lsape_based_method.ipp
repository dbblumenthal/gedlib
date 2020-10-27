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
 * @file  lsape_based_method.ipp
 * @brief ged::LSAPEBasedMethod class defeinition.
 */

#ifndef SRC_METHODS_LSAPE_BASED_METHOD_IPP_
#define SRC_METHODS_LSAPE_BASED_METHOD_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
~LSAPEBasedMethod() {}

template<class UserNodeLabel, class UserEdgeLabel>
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
LSAPEBasedMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
lsape_model_{LSAPESolver::ECBP},
greedy_method_{LSAPESolver::BASIC},
enumeration_method_{LSAPESolver::DISSIMILAR},
compute_lower_bound_{true},
solve_optimally_{true},
num_threads_{1},
centrality_method_{NONE},
centrality_weight_{0.7},
centralities_(),
max_num_solutions_{1} {}

// === Definition of public member functions that extend the public interface provided by GEDMethod.

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
populate_instance(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance) {
	lsape_instance.resize(g.num_nodes() + 1, h.num_nodes() + 1);
	if (not this->initialized_) {
		lsape_pre_graph_init_(true);
		init_graph_(g);
		init_graph_(h);
		lsape_default_post_graph_init_();
	}
	lsape_populate_instance_(g, h, lsape_instance);
	lsape_instance(g.num_nodes(), h.num_nodes()) = 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
populate_instance_and_run_as_util(const GEDGraph & g, const GEDGraph & h, Result & result, DMatrix & lsape_instance) {

	// Populate the LSAPE instance and set up the solver.
	populate_instance(g, h, lsape_instance);
	LSAPESolver lsape_solver(&lsape_instance);

	// Solve the LSAPE instance.
	if (solve_optimally_) {
		lsape_solver.set_model(lsape_model_);
        lsape_solver.set_enumeration_method(enumeration_method_);
	}
	else {
		lsape_solver.set_greedy_method(greedy_method_);
	}
	lsape_solver.solve(max_num_solutions_);

	// Compute and store lower and upper bound
	if (compute_lower_bound_ and solve_optimally_) {
		result.set_lower_bound(lsape_solver.minimal_cost() * lsape_lower_bound_scaling_factor_(g, h));
	}

	for (std::size_t solution_id{0}; solution_id < lsape_solver.num_solutions(); solution_id++) {
		std::size_t index_node_map{result.add_node_map(g.num_nodes(), h.num_nodes())};
		util::construct_node_map_from_solver(lsape_solver, result.node_map(index_node_map), solution_id);
		this->ged_data_.compute_induced_cost(g, h, result.node_map(index_node_map));
	}

	// Add centralities and reoptimize
	if (centrality_weight_ > 0.0 and centrality_method_ != NONE) {
		add_centralities_(g, h, lsape_instance);
		lsape_solver.solve(max_num_solutions_);
		for (std::size_t solution_id{0}; solution_id < lsape_solver.num_solutions(); solution_id++) {
			std::size_t index_node_map{result.add_node_map(g.num_nodes(), h.num_nodes())};
			util::construct_node_map_from_solver(lsape_solver, result.node_map(index_node_map), solution_id);
			this->ged_data_.compute_induced_cost(g, h, result.node_map(index_node_map));
		}
	}

	// Sort the node maps an set the upper bound.
	result.sort_node_maps_and_set_upper_bound(max_num_solutions_);
}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
	lsape_pre_graph_init_(false);
	for (auto graph = this->ged_data_.begin(); graph != this->ged_data_.end(); graph++) {
		init_graph_(*graph);
	}
	lsape_init_();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {
	DMatrix lsape_instance;
	populate_instance_and_run_as_util(g, h, result, lsape_instance);
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	bool is_valid_option{false};
	if (option == "threads") {
		try {
			num_threads_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		if (num_threads_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "lsape-model") {
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
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option lsape-model. Usage: options = \"[--lsape-model ECBP|EBP|FLWC|FLCC|FBP|SFBP|FBP0] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "greedy-method") {
		if (arg == "BASIC") {
			greedy_method_ = LSAPESolver::BASIC;
		}
		else if (arg == "REFINED") {
			greedy_method_ = LSAPESolver::REFINED;
		}
		else if (arg == "LOSS") {
			greedy_method_ = LSAPESolver::LOSS;
		}
		else if (arg == "BASIC_SORT") {
			greedy_method_ = LSAPESolver::BASIC_SORT;
		}
		else if (arg == "INT_BASIC_SORT") {
			greedy_method_ = LSAPESolver::INT_BASIC_SORT;
		}
		else {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option greedy-method. Usage: options = \"[--greedy-method BASIC|REFINED|LOSS|BASIC_SORT|INT_BASIC_SORT] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "enumeration-method") {
	    if (arg == "BASELINE") {
	        enumeration_method_ = LSAPESolver::BASELINE;
	    }
	    else if (arg == "DISSIMILAR") {
            enumeration_method_ = LSAPESolver::DISSIMILAR;
	    }
	    else {
            throw ged::Error(std::string("Invalid argument ") + arg  + " for option enumeration-method. Usage: options = \"[--enumeration-method BASELINE|DISSIMILAR] [...]\"");
	    }
        is_valid_option = true;
	}
	else if (option == "optimal") {
		if (arg == "TRUE") {
			solve_optimally_ = true;
		}
		else if (arg == "FALSE") {
			solve_optimally_ = false;
		}
		else {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option optimal. Usage: options = \"[--optimal TRUE|FALSE] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "centrality-method") {
		if (arg == "NONE") {
			centrality_method_ = NONE;
		}
		else if (arg == "DEGREE") {
			centrality_method_ = DEGREE;
		}
		else if (arg == "EIGENVECTOR") {
			centrality_method_ = EIGENVECTOR;
		}
		else if (arg == "PAGERANK") {
			centrality_method_ = PAGERANK;
		}
		else {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option centrality-method. Usage: options = \"[--centrality-method NONE|DEGREE|EIGENVECTOR|PAGERANK] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "centrality-weight") {
		try {
			centrality_weight_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument ") + arg + " for option centrality-weight. Usage: options = \"[--centrality-weight <convertible to double between 0 and 1>] [...]");
		}
		if (centrality_weight_ < 0.0 or centrality_weight_ > 1.0) {
			throw Error(std::string("Invalid argument ") + arg + " for option centrality-weight. Usage: options = \"[--centrality-weight <convertible to double between 0 and 1>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "max-num-solutions") {
		if (arg == "ALL") {
			max_num_solutions_ = -1;
		}
		else {
			try {
				max_num_solutions_ = std::stoi(arg);
			}
			catch (...) {
				throw Error(std::string("Invalid argument ") + arg + " for option  max-num-solutions. Usage: options = \"[--max-num-solutions ALL|<convertible to int greater 0>] [...]");
			}
			if (max_num_solutions_ < 1) {
				throw Error(std::string("Invalid argument ") + arg + " for option  max-num-solutions. Usage: options = \"[--max-num-solutions ALL|<convertible to int greater 0>] [...]");
			}
		}
		is_valid_option = true;
	}
	is_valid_option = is_valid_option or lsape_parse_option_(option, arg);
	return is_valid_option;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	if (lsape_valid_options_string_() == "") {
		return "[--lsape-model <arg>] [--threads <arg>] [--greedy-method <arg>] [--enumeration-method <arg>] [--optimal <arg>] [--centrality-method <arg>] [--centrality-weight <arg>] [--max-num-solutions <arg>]";
	}
	return lsape_valid_options_string_() + " [--lsape-model <arg>] [--greedy-method <arg>] [--enumeration-method <arg>] [--optimal <arg>] [--centrality-method <arg>] [--centrality-weight <arg>] [--max-num-solutions <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	lsape_model_ = LSAPESolver::ECBP;
	greedy_method_ = LSAPESolver::BASIC;
	enumeration_method_ = LSAPESolver::DISSIMILAR;
	solve_optimally_ = true;
	centrality_method_ = NONE;
	centrality_weight_ = 0.7;
	max_num_solutions_ = 1;
	num_threads_ = 1;
	lsape_set_default_options_();
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
init_graph_(const GEDGraph & graph) {
	if (centrality_method_ != NONE) {
		init_centralities_(graph);
	}
	lsape_init_graph_(graph);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
init_centralities_(const GEDGraph & graph) {
	centralities_[graph.id()] = std::vector<double>(graph.num_nodes(), 0.0);
	if (centrality_method_ == DEGREE) {
		for (std::size_t row{0}; row < graph.num_nodes(); row++) {
			centralities_.at(graph.id())[row] = static_cast<double>(graph.degree(row));
		}
		return;
	}
	DMatrix adj_matrix(graph.num_nodes(), graph.num_nodes());
	util::init_adj_matrix(graph, adj_matrix);
	double eigenvalue;
	std::vector<double> eigenvector;
	compute_eigenvector_with_largest_eigenvalue_(adj_matrix, eigenvector, eigenvalue);
	for (std::size_t row{0}; row < adj_matrix.num_rows(); row++) {
		for (std::size_t col{0}; col < adj_matrix.num_cols(); col++) {
			double summand{adj_matrix(row, col) * eigenvector.at(col)};
			if (centrality_method_ == PAGERANK) {
				summand /= std::max(1.0, static_cast<double>(graph.degree(col)));
				summand += 1.0;
			}
			centralities_.at(graph.id())[row] += summand;
		}
		if (centrality_method_ == PAGERANK) {
			centralities_.at(graph.id())[row] *= 0.85;
		}
		else {
			centralities_.at(graph.id())[row] /= eigenvalue;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
add_centralities_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {
	// substitution
	master_problem *= (1 - centrality_weight_);
#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row = 0; row < master_problem.num_rows(); row++) {
		for (std::size_t col = 0; col < master_problem.num_cols(); col++) {
			if (row < g.num_nodes() and col < h.num_nodes()) {
				master_problem(row, col) += centrality_weight_ * std::fabs(static_cast<double>(centralities_.at(g.id()).at(row) - centralities_.at(h.id()).at(col)));
			}
			else if (row < g.num_nodes()) {
				master_problem(row, h.num_nodes()) += centrality_weight_ * centralities_.at(g.id()).at(row);
			}
			else if (col < h.num_nodes()) {
				master_problem(g.num_nodes(), col) += centrality_weight_ * centralities_.at(h.id()).at(col);
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
compute_eigenvector_with_largest_eigenvalue_(const DMatrix & symmetric_matrix, std::vector<double> & eigenvector, double & eigenvalue) {
	if (symmetric_matrix.num_rows() != symmetric_matrix.num_cols()) {
		return;
	}
	Eigen::EigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> eigen_solver(symmetric_matrix.matrix());
	Eigen::VectorXd eigen_values = eigen_solver.eigenvalues().real();
	std::size_t pos;
	eigenvalue = eigen_values.maxCoeff(&pos);
	Eigen::VectorXd eigen_vector = eigen_solver.eigenvectors().col(pos).real();
	eigenvector.resize(eigen_vector.rows());
	for (std::size_t row{0}; row < static_cast<std::size_t>(eigen_vector.rows()); row++) {
		eigenvector[row] = eigen_vector(row);
	}
}

// === Default definitions of private virtual member functions to be overridden by derived classes. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_init_() {}

template<class UserNodeLabel, class UserEdgeLabel>
bool
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_parse_option_(const std::string & option, const std::string & arg) {
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	return "";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_pre_graph_init_(bool called_at_runtime) {}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_default_post_graph_init_() {}

template<class UserNodeLabel, class UserEdgeLabel>
double
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>::
lsape_lower_bound_scaling_factor_(const GEDGraph & g, const GEDGraph & h) {
	return 1.0;
}

}

#endif /* SRC_METHODS_LSAPE_BASED_METHOD_IPP_ */
