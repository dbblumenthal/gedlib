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
 * @file  branch_tight.ipp
 * @brief ged::BranchTight class definition.
 */

#ifndef SRC_METHODS_BRANCH_TIGHT_IPP_
#define SRC_METHODS_BRANCH_TIGHT_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchTight<UserNodeLabel, UserEdgeLabel>::
~BranchTight() {}

template<class UserNodeLabel, class UserEdgeLabel>
BranchTight<UserNodeLabel, UserEdgeLabel>::
BranchTight(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
max_itrs_{20},
range_{0.0},
epsilon_{0.01},
upper_bound_option_{BEST},
naive_regularization_{true},
num_threads_{1},
time_limit_in_sec_{0.0} {}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	max_itrs_ = 20;
	range_ = 0.0;
	epsilon_ = 0.01;
	upper_bound_option_ = BEST;
	naive_regularization_ = true;
	num_threads_= 1;
	time_limit_in_sec_ = 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
BranchTight<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	return "[--iterations <arg>] [--range <arg>] [--epsilon <arg>] [--upper-bound <arg>] [--regularize <arg>] [--threads <arg>] [--time-limit <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchTight<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "iterations") {
		try {
			max_itrs_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option iterations. Usage: options = \"[--iterations <convertible to int>] [...]");
		}
		if ((max_itrs_ <= 0) and (epsilon_ <= 0) and (time_limit_in_sec_ <= 0)) {
			throw Error("Switching off all termination criteria that guarantee termination is illegal.");
		}
		return true;
	}
	else if (option == "time-limit") {
		try {
			time_limit_in_sec_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option time-limit. Usage: options = \"[--time-limit <convertible to double>] [...]");
		}
		if ((max_itrs_ <= 0) and (epsilon_ <= 0) and (time_limit_in_sec_ <= 0)) {
			throw Error("Switching off all termination criteria that guarantee termination is illegal.");
		}
		return true;
	}
	else if (option == "range") {
		try {
			range_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option range. Usage: options = \"[--range <convertible to double>] [...]");
		}
		return true;
	}
	else if (option == "epsilon") {
		try {
			epsilon_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option epsilon. Usage: options = \"[--epsilon <convertible to double>] [...]");
		}
		if ((max_itrs_ <= 0) and (epsilon_ <= 0) and (time_limit_in_sec_ <= 0)) {
			throw Error("Switching off all termination criteria that guarantee termination is illegal.");
		}
		return true;
	}
	else if (option == "upper-bound") {
		if (arg == "FIRST") {
			upper_bound_option_ = FIRST;
		}
		else if (arg == "LAST") {
			upper_bound_option_ = LAST;
		}
		else if (arg == "BEST") {
			upper_bound_option_ = BEST;
		}
		else if (arg == "NO") {
			upper_bound_option_ = NO;
		}
		else {
			throw ged::Error(std::string("Invalid argument ") + arg  + " for option upper-bound. Usage: options = \"[--upper-bound NO|FIRST|LAST|BEST]\"");
		}
		return true;
	}
	else if (option == "regularize") {
		if (arg == "K-FACTOR") {
			naive_regularization_ = false;
		}
		else if (arg == "NAIVE") {
			naive_regularization_ = true;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option range. Usage: options = \"[--regularize NAIVE|K-FACTOR] [...]");
		}
		return true;
	}
	else if (option == "threads") {
		try {
			num_threads_ = std::stoi(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		if (max_itrs_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		return true;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & gg, const GEDGraph & hh, Result & result) {

	Timer timer(time_limit_in_sec_);

	// Regularize the graphs g and h.
	GEDGraph g(gg);
	GEDGraph h(hh);
	std::size_t degree{};
	if (naive_regularization_) {
		degree = naive_regularize_(g, h);
	}
	else {
		degree = regularize_(g, h);
	}


	// Initialize the subproblems, the node costs, and the weights for updating the subproblems.
	SubproblemSolvers_ subproblem_solvers(static_cast<std::size_t>(g.num_nodes()), degree);
	init_subproblems_(g, h, subproblem_solvers);
	DMatrix node_costs(static_cast<std::size_t>(g.num_nodes()), static_cast<std::size_t>(g.num_nodes()));
	init_node_costs_(g, h, node_costs);
	Weights_ weights(static_cast<std::size_t>(g.num_nodes()));

	// The main loop.
	DMatrix master_problem(g.num_nodes(), g.num_nodes());
	LSAPSolver master_problem_solver(&master_problem);
	double last_improvement{std::numeric_limits<double>::max()};
	for (std::size_t current_itr{1}; not termination_criterion_met_(current_itr, last_improvement, result); current_itr++) {

		if ((current_itr > 1) and timer.expired()) {
			break;
		}

		// Update and solve the subproblems in parallel.
		if (current_itr > 1) {
			update_weights_(master_problem_solver, degree, subproblem_solvers, weights);
			update_subproblem_costs_(weights, degree, subproblem_solvers);
		}
		subproblem_solvers.solve(num_threads_);

		// Update and solve the master problem.
		update_master_problem_costs_(subproblem_solvers, node_costs, master_problem);
		master_problem_solver.solve();

		// Update the lower bound bound, the improvement ratio, and the upper bound, if necessary.
		if (current_itr > 1) {
			if (result.lower_bound() < 0.00001) {
				break;
			}
			last_improvement = (master_problem_solver.minimal_cost() - result.lower_bound()) / result.lower_bound();
		}
		result.set_lower_bound(master_problem_solver.minimal_cost());
		if (upper_bound_option_ == BEST or (upper_bound_option_ == FIRST and current_itr == 1)) {
			std::size_t index_node_map{result.add_node_map(g.num_nodes(), h.num_nodes())};
			util::construct_node_map_from_solver(master_problem_solver, result.node_map(index_node_map));
			if (result.is_non_redundant_node_map(index_node_map)) {
				this->ged_data_.compute_induced_cost(g, h, result.node_map(index_node_map));
				result.sort_node_maps_and_set_upper_bound();
			}
		}
	}

	// Store the last upper bound if the option upper-bound is set to LAST.
	if (upper_bound_option_ == LAST) {
		std::size_t index_node_map{result.add_node_map(g.num_nodes(), h.num_nodes())};
		util::construct_node_map_from_solver(master_problem_solver, result.node_map(index_node_map));
		this->ged_data_.compute_induced_cost(g, h, result.node_map(index_node_map));
		result.sort_node_maps_and_set_upper_bound();
	}
}

// === Definitions of private helper functions. ===

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchTight<UserNodeLabel, UserEdgeLabel>::
termination_criterion_met_(const std::size_t & current_itr, const double & last_improvement, Result & result) {
	if ((result.upper_bound() >= 0) and (result.lower_bound() >= result.upper_bound())) {
		return true;
	}
	if ((max_itrs_ > 0) and (max_itrs_ < current_itr)) {
		return true;
	}
	if ((epsilon_ > 0) and (epsilon_ > last_improvement)) {
		return true;
	}
	if ((range_ > 0) and (range_ < result.lower_bound())) {
		return true;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
init_subproblems_(const GEDGraph & g, const GEDGraph & h, SubproblemSolvers_ & subproblems_solver) const {
	for (std::size_t row_master{0}; row_master < subproblems_solver.get_size(); row_master++) {

		// Collect the edges that are incident to the node that corresponds to row row_master in the master problem.
		auto incident_edges_i = g.incident_edges(row_master);

		for (std::size_t col_master{0}; col_master < subproblems_solver.get_size(); col_master++) {

			// Collect the edges that are incident to the node that corresponds to column col_master in the master problem.
			auto incident_edges_k = h.incident_edges(col_master);

			// Initialize row and collumn indices.
			std::size_t row_sub{0};
			for (auto ij = incident_edges_i.first; ij != incident_edges_i.second; ij++) {
				subproblems_solver.rows_subproblem_to_master(row_master, col_master)[row_sub++] = g.head(*ij);
			}
			std::size_t col_sub{0};
			for (auto kl = incident_edges_k.first; kl != incident_edges_k.second; kl++) {
				subproblems_solver.cols_subproblem_to_master(row_master, col_master)[col_sub++] = h.head(*kl);
			}

			// Collect edge relabelling costs.
			row_sub = 0;
			for (auto ij = incident_edges_i.first; ij != incident_edges_i.second; ij++) {
				col_sub = 0;
				for (auto kl = incident_edges_k.first; kl != incident_edges_k.second; kl++) {
					subproblems_solver.subproblem(row_master, col_master)(row_sub, col_sub++) = this->ged_data_.edge_cost(g.get_edge_label(*ij), h.get_edge_label(*kl)) * 0.5;
				}
				row_sub++;
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
update_master_problem_costs_(const SubproblemSolvers_ & subproblems_solver, const DMatrix & node_costs, DMatrix & master_problem) const {
	for (std::size_t row_master{0}; row_master < subproblems_solver.get_size(); row_master++) {
		for (std::size_t col_master{0}; col_master < subproblems_solver.get_size(); col_master++) {
			master_problem(row_master, col_master) = node_costs(row_master, col_master) + subproblems_solver.solver(row_master, col_master).minimal_cost();
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
init_node_costs_(const GEDGraph & g, const GEDGraph & h, DMatrix & node_costs) const {
	for (std::size_t row_master{0}; row_master < node_costs.num_rows(); row_master++) {
		for (std::size_t col_master{0}; col_master < node_costs.num_cols(); col_master++) {
			node_costs(row_master, col_master) = this->ged_data_.node_cost(g.get_node_label(row_master), h.get_node_label(col_master));
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
update_weights_(const LSAPSolver & master_problem_solver, std::size_t degree, const SubproblemSolvers_ & subproblems_solver, Weights_ & weights) const {
	double delta{1 / static_cast<double>(degree)};
	double weight{0.0};
	for (std::size_t row_master{0}; row_master < subproblems_solver.get_size(); row_master++) {
		for (std::size_t col_master{0}; col_master < subproblems_solver.get_size(); col_master++) {
			for (std::size_t row_sub{0}; row_sub < degree; row_sub++) {
				std::size_t row_sub_in_master{subproblems_solver.row_in_master(row_master, col_master, row_sub)};
				for (std::size_t col_sub{0}; col_sub < degree; col_sub++) {
					std::size_t col_sub_in_master{subproblems_solver.col_in_master(row_master, col_master, col_sub)};
					weight = subproblems_solver.solver(row_master, col_master).get_slack(row_sub, col_sub) + master_problem_solver.get_slack(row_master, col_master) * delta;
					weights.set_weight(row_master, col_master, row_sub_in_master, col_sub_in_master, weight);
				}
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
update_subproblem_costs_(const Weights_ & weights, std::size_t degree, SubproblemSolvers_ & subproblems_solver) const {
	double added_weight{0.0};
	for (std::size_t row_master{0}; row_master < subproblems_solver.get_size(); row_master++) {
		for (std::size_t col_master{0}; col_master < subproblems_solver.get_size(); col_master++) {
			for (std::size_t row_sub{0}; row_sub < degree; row_sub++) {
				std::size_t row_sub_in_master{subproblems_solver.row_in_master(row_master, col_master, row_sub)};
				for (std::size_t col_sub{0}; col_sub < degree; col_sub++) {
					std::size_t col_sub_in_master{subproblems_solver.col_in_master(row_master, col_master, col_sub)};
					added_weight = weights.get_weight(row_sub_in_master, col_sub_in_master, row_master, col_master);
					added_weight -= weights.get_weight(row_master, col_master, row_sub_in_master, col_sub_in_master);
					subproblems_solver.subproblem(row_master, col_master)(row_sub, col_sub) += added_weight;
				}
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BranchTight<UserNodeLabel, UserEdgeLabel>::
regularize_(GEDGraph & g, GEDGraph & h) const {

	// Add dummy nodes to smaller graph or to both graphs if the costs are not quasimetric.
	if (this->ged_data_.quasimetric_costs(g, h)) {
		fill_up_smaller_graph_(g, h);
	}
	else {
		fill_up_both_graphs_(g, h);
	}

	// Build complement graphs.
	GEDGraph complement_g{};
	GEDGraph::NodeNodeMap complement_to_g;
	construct_complement_graph_(g, complement_g, complement_to_g);
	GEDGraph complement_h{};
	GEDGraph::NodeNodeMap complement_to_h;
	construct_complement_graph_(h, complement_h, complement_to_h);

	// Initialize k-factors.
	GEDGraph::NodeSizeTMap complement_graph_g_to_k_factor, g_to_k_factor;
	std::size_t size_k_factor_g{0};
	for (auto complement_node = complement_g.nodes().first; complement_node != complement_g.nodes().second; complement_node++) {
		complement_graph_g_to_k_factor[*complement_node] = size_k_factor_g;
		g_to_k_factor[complement_to_g.at(*complement_node)] = size_k_factor_g++;
	}
	KFactor_ k_factor_g(size_k_factor_g);
	GEDGraph::NodeSizeTMap complement_graph_h_to_k_factor, h_to_k_factor;
	std::size_t size_k_factor_h{0};
	for (auto complement_node = complement_h.nodes().first; complement_node != complement_h.nodes().second; complement_node++) {
		complement_graph_h_to_k_factor[*complement_node] = size_k_factor_h;
		h_to_k_factor[complement_to_h.at(*complement_node)] = size_k_factor_h++;
	}
	KFactor_ k_factor_h(size_k_factor_h);

	// Find maximum k s.t. both complement graphs have a k-factor.
	std::size_t num_nodes(g.num_nodes());
	std::size_t max_k(num_nodes - 1 - std::max(g.maxdeg(), h.maxdeg()));
	std::size_t k{max_k};
	for (; k >= 1; k--) {
		if (not compute_k_factor_(complement_g, k, complement_graph_g_to_k_factor, k_factor_g)) {
			k_factor_g.clear_edges();
			continue;
		}
		if (not compute_k_factor_(complement_h, k, complement_graph_h_to_k_factor, k_factor_h)) {
			k_factor_g.clear_edges();
			k_factor_h.clear_edges();
			continue;
		}
		break;
	}

	// Regularize g and h by adding those edges that are in the complements but not in the k-factors.
	regularize_from_k_factor_(k_factor_g, g_to_k_factor, g);
	regularize_from_k_factor_(k_factor_h, h_to_k_factor, h);

	// Return the degree of all nodes in the regularized graphs g and h.
	return (num_nodes - 1 - k);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BranchTight<UserNodeLabel, UserEdgeLabel>::
naive_regularize_(GEDGraph & g, GEDGraph & h) const {

	// Add dummy nodes to smaller graph or to both graphs if the costs are not quasimetric.
	if (this->ged_data_.quasimetric_costs()) {
		fill_up_smaller_graph_(g, h);
	}
	else {
		fill_up_both_graphs_(g, h);
	}

	// Add missing edges to g.
	for (auto node_1 = g.nodes().first; node_1 != g.nodes().second; node_1++) {
		for (auto node_2 = node_1 + 1; node_2 != g.nodes().second; node_2++) {
			if (not g.is_edge(*node_1, *node_2)) {
				GEDGraph::EdgeID new_edge{g.add_edge(*node_1, *node_2)};
				g.set_label(new_edge, dummy_label());
			}
		}
	}
	g.setup_adjacency_matrix();

	// add missing edges to h.
	for (auto node_1 = h.nodes().first; node_1 != h.nodes().second; node_1++) {
		for (auto node_2 = node_1 + 1; node_2 != h.nodes().second; node_2++) {
			if (not h.is_edge(*node_1, *node_2)) {
				GEDGraph::EdgeID new_edge{h.add_edge(*node_1, *node_2)};
				h.set_label(new_edge, dummy_label());
			}
		}
	}
	h.setup_adjacency_matrix();

	return (g.num_nodes() - 1);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
regularize_from_k_factor_(const KFactor_ & k_factor, const GEDGraph::NodeSizeTMap & graph_to_k_factor, GEDGraph & graph) const {
	for (auto node_1 = graph.nodes().first; node_1 != graph.nodes().second; node_1++) {
		for (auto node_2 = node_1 + 1; node_2 != graph.nodes().second; node_2++) {
			if (graph.is_edge(*node_1, *node_2)) {
				continue;
			}
			if (not k_factor.contains_edge(graph_to_k_factor.at(*node_1), graph_to_k_factor.at(*node_2))) {
				GEDGraph::EdgeID new_edge{graph.add_edge(*node_1, *node_2)};
				graph.set_label(new_edge, dummy_label());
			}
		}
	}
	graph.setup_adjacency_matrix();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
fill_up_smaller_graph_(GEDGraph & g, GEDGraph & h) const {
	if (g.num_nodes() < h.num_nodes()) {
		for (std::size_t counter{g.num_nodes()}; counter < h.num_nodes(); counter++) {
			GEDGraph::NodeID new_node{g.add_node()};
			g.set_label(new_node, dummy_label());
		}
		g.setup_adjacency_matrix();
	}
	else if (g.num_nodes() > h.num_nodes()) {
		for (std::size_t counter{h.num_nodes()}; counter < g.num_nodes(); counter++) {
			GEDGraph::NodeID new_node{h.add_node()};
			h.set_label(new_node, dummy_label());
		}
		h.setup_adjacency_matrix();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
fill_up_both_graphs_(GEDGraph & g, GEDGraph & h) const {
	std::size_t g_num_nodes{g.num_nodes()};
	for (std::size_t counter{0}; counter < h.num_nodes(); counter++) {
		GEDGraph::NodeID new_node{g.add_node()};
		g.set_label(new_node, dummy_label());
	}
	g.setup_adjacency_matrix();
	for (std::size_t counter{0}; counter < g_num_nodes; counter++) {
		GEDGraph::NodeID new_node{h.add_node()};
		h.set_label(new_node, dummy_label());
	}
	h.setup_adjacency_matrix();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
construct_complement_graph_(const GEDGraph & graph, GEDGraph & complement_graph, GEDGraph::NodeNodeMap & complement_to_graph) const {
	GEDGraph::NodeNodeMap graph_to_complement;
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		GEDGraph::NodeID complement_node{complement_graph.add_node()};
		graph_to_complement[*node] = complement_node;
		complement_to_graph[complement_node] = *node;
	}
	for (auto node_1 = graph.nodes().first; node_1 != graph.nodes().second; node_1++) {
		for (auto node_2 = node_1 + 1; node_2 != graph.nodes().second; node_2++) {
			if (not graph.is_edge(*node_1, *node_2)) {
				complement_graph.add_edge(graph_to_complement.at(*node_1),graph_to_complement.at(*node_2));
			}
		}
	}
	complement_graph.setup_adjacency_matrix();
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchTight<UserNodeLabel, UserEdgeLabel>::
compute_k_factor_(const GEDGraph & complement_graph, std::size_t k, const GEDGraph::NodeSizeTMap & complement_graph_to_k_factor, KFactor_ & k_factor) const {

	// Initialize output variables for construct_transformed_complement_graph_().
	GEDGraph transformed_complement_graph{};
	GEDGraph::NodeNodeMap transformed_to_original_nodes;

	/*
	 * Call construct_transformed_complement_graph_() and return false if the
	 * construction failed because k is larger than the minimum degree in the
	 * complement graph.
	 */
	if (not construct_transformed_complement_graph_(complement_graph, k, transformed_complement_graph, transformed_to_original_nodes)) {
		return false;
	}

	// Initialize matching.
	GEDGraph::NodeNodeMap matching;

	/*
	 * Use the boost implementation of edmonds blossom algorithm to compute a
	 * maximum matching in transformed complement graph. Reconstruct k-factor
	 * for the complement graph and return true if the matching is perfect.
	 * Otherwise, return false.
	 */
	boost::edmonds_maximum_cardinality_matching(transformed_complement_graph.get_adjacency_list(), MateMap_(matching));
	for (auto assignment = matching.begin(); assignment != matching.end(); assignment++) {
		if (assignment->second == GEDGraph::dummy_node()) {
			return false;
		}
		auto original_node = transformed_to_original_nodes.find(assignment->first);
		auto assigned_original_node = transformed_to_original_nodes.find(assignment->second);
		if (original_node != transformed_to_original_nodes.end() and assigned_original_node != transformed_to_original_nodes.end()) {
			k_factor.add_edge(complement_graph_to_k_factor.at(original_node->second), complement_graph_to_k_factor.at(assigned_original_node->second));
		}
	}
	return true;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchTight<UserNodeLabel, UserEdgeLabel>::
construct_transformed_complement_graph_(const GEDGraph & complement_graph, std::size_t k, GEDGraph & transformed_complement_graph, GEDGraph::NodeNodeMap & transformed_to_original_nodes) const {

	// colect the degrees of the transformed graph's nodes and return false is mindeg < k
	GEDGraph::NodeSizeTMap degrees;
	for (auto node = complement_graph.nodes().first; node != complement_graph.nodes().second; node++) {
		degrees[*node] = complement_graph.degree(*node);
		if (degrees[*node] < k) {
			return false;
		}
	}

	// check if to use type 1 gadgets or type 2 gadgets
	std::map<GEDGraph::NodeID,bool> type_1_gadget;
	for (auto node_deg_pair = degrees.begin(); node_deg_pair != degrees.end(); node_deg_pair++) {
		type_1_gadget[node_deg_pair->first] = (k < (node_deg_pair->second / 2));
	}

	// the IDs of the core, inner, and outer nodes
	std::map<GEDGraph::NodeID, std::vector<GEDGraph::NodeID>> core_nodes;
	std::map<GEDGraph::NodeID, std::vector<GEDGraph::NodeID>> inner_nodes;
	std::map<GEDGraph::NodeID, std::vector<GEDGraph::NodeID>> outer_nodes;
	for (auto node = complement_graph.nodes().first; node != complement_graph.nodes().second; node++) {
		core_nodes[*node] = std::vector<GEDGraph::NodeID>();
		inner_nodes[*node] = std::vector<GEDGraph::NodeID>();
		outer_nodes[*node] = std::vector<GEDGraph::NodeID>();
	}

	// add nodes to the transformed graph
	for (auto node = complement_graph.nodes().first; node != complement_graph.nodes().second; node++) {
		if (type_1_gadget.at(*node)) {
			for (std::size_t counter{0}; counter < k; counter++) {
				core_nodes[*node].push_back(transformed_complement_graph.add_node());
			}
			for (std::size_t counter{0}; counter < degrees.at(*node); counter++) {
				inner_nodes[*node].push_back(transformed_complement_graph.add_node());
			}
		}
		else {
			for (std::size_t counter{0}; counter < degrees.at(*node) - k; counter++) {
				core_nodes[*node].push_back(transformed_complement_graph.add_node());
			}
		}
		for (std::size_t counter{0}; counter < degrees.at(*node); counter++) {
			GEDGraph::NodeID transformed_node{transformed_complement_graph.add_node()};
			outer_nodes[*node].push_back(transformed_node);
			transformed_to_original_nodes[transformed_node] = *node;
		}
	}

	// add core and inner edges to the transformed graph
	for (auto node = complement_graph.nodes().first; node != complement_graph.nodes().second; node++) {
		if (type_1_gadget.at(*node)) {
			for (auto inner_node = inner_nodes.at(*node).begin(), outer_node = outer_nodes.at(*node).begin(); inner_node != inner_nodes.at(*node).end(); inner_node++, outer_node++) {
				transformed_complement_graph.add_edge(*inner_node, *outer_node);
				for (auto core_node = core_nodes.at(*node).begin(); core_node != core_nodes.at(*node).end(); core_node++) {
					transformed_complement_graph.add_edge(*core_node, *inner_node);
				}
			}
		}
		else {
			for (auto core_node = core_nodes.at(*node).begin(); core_node != core_nodes.at(*node).end(); core_node++) {
				for (auto outer_node = outer_nodes.at(*node).begin(); outer_node != outer_nodes.at(*node).end(); outer_node++) {
					transformed_complement_graph.add_edge(*core_node, *outer_node);
				}
			}
		}
	}

	// add outer edges to the transformed graph
	std::map<GEDGraph::NodeID, std::vector<GEDGraph::NodeID>::iterator> outer_nodes_itr;
	for (auto node = complement_graph.nodes().first; node != complement_graph.nodes().second; node++) {
		outer_nodes_itr[*node] = outer_nodes.at(*node).begin();
	}
	for (auto node_1 = complement_graph.nodes().first; node_1 != complement_graph.nodes().second; node_1++) {
		for (auto node_2 = node_1 +1; node_2 != complement_graph.nodes().second; node_2++) {
			if (complement_graph.is_edge(*node_1, *node_2)) {
				transformed_complement_graph.add_edge(*outer_nodes_itr.at(*node_1)++, *outer_nodes_itr.at(*node_2)++);
			}
		}
	}
	return true;
}

// === Definition of private class KFactor_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchTight<UserNodeLabel, UserEdgeLabel>::
KFactor_ ::
KFactor_(LabelID num_nodes) :
num_nodes_{num_nodes},
is_edge_(num_nodes_ * num_nodes_, false) {}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchTight<UserNodeLabel, UserEdgeLabel>::
KFactor_ ::
contains_edge(std::size_t id_i, std::size_t id_k) const {
	return is_edge_.at(id_i + id_k * num_nodes_);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BranchTight<UserNodeLabel, UserEdgeLabel>::
KFactor_ ::
num_nodes() const {
	return num_nodes_;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
KFactor_ ::
add_edge(std::size_t id_i, std::size_t id_k) {
	is_edge_[id_i + id_k * num_nodes_] = true;
	is_edge_[id_k + id_i * num_nodes_] = true;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
KFactor_ ::
clear_edges() {
	for (auto edge_info = is_edge_.begin(); edge_info != is_edge_.end(); edge_info++) {
		*edge_info = false;
	}
}

// === Definition of private class Weights_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchTight<UserNodeLabel, UserEdgeLabel>::
Weights_ ::
Weights_(std::size_t size_master) :
size_master_{size_master},
weights_(size_master_ * size_master_ * size_master_ * size_master_, 0.0){}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchTight<UserNodeLabel, UserEdgeLabel>::
Weights_ ::
get_weight(std::size_t row_in_master, std::size_t col_in_master, std::size_t row_sub_in_master, std::size_t col_sub_in_master) const {
	return weights_.at(row_in_master + size_master_ * col_in_master + size_master_ * size_master_ * row_sub_in_master + size_master_ * size_master_ * size_master_ * col_sub_in_master);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
Weights_ ::
set_weight(std::size_t row_in_master, std::size_t col_in_master, std::size_t row_sub_in_master, std::size_t col_sub_in_master, double weight) {
	weights_[row_in_master + size_master_ * col_in_master + size_master_ * size_master_ * row_sub_in_master + size_master_ * size_master_ * size_master_ * col_sub_in_master] = weight;
}

// === Definition of private class SubproblemSolvers_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
SubproblemSolvers_(std::size_t size_master, std::size_t degree) :
size_master_{size_master},
degree_{degree},
subproblems_(size_master_ * size_master_, DMatrix(degree_, degree_)),
subproblem_solvers_(),
rows_sub_to_master_(size_master_ * size_master_),
cols_sub_to_master_(size_master_ * size_master_) {
	for (const auto & subproblem : subproblems_) {
		subproblem_solvers_.emplace_back(&subproblem);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
LSAPSolver &
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
solver(std::size_t row_in_master, std::size_t col_in_master) {
	return subproblem_solvers_.at(row_in_master + size_master_ * col_in_master);
}

template<class UserNodeLabel, class UserEdgeLabel>
const LSAPSolver &
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
solver(std::size_t row_in_master, std::size_t col_in_master) const {
	return subproblem_solvers_.at(row_in_master + size_master_ * col_in_master);
}

template<class UserNodeLabel, class UserEdgeLabel>
DMatrix &
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
subproblem(std::size_t row_in_master, std::size_t col_in_master) {
	return subproblems_.at(row_in_master + size_master_ * col_in_master);
}

template<class UserNodeLabel, class UserEdgeLabel>
const DMatrix &
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
subproblem(std::size_t row_in_master, std::size_t col_in_master) const {
	return subproblems_.at(row_in_master + size_master_ * col_in_master);
}

template<class UserNodeLabel, class UserEdgeLabel>
typename BranchTight<UserNodeLabel, UserEdgeLabel>::SizeTSizeTMap_ &
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
rows_subproblem_to_master(std::size_t row_in_master, std::size_t col_in_master) {
	return rows_sub_to_master_.at(row_in_master + size_master_ * col_in_master);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
row_in_master(std::size_t row_in_master, std::size_t col_in_master, std::size_t row_in_subproblem) const {
	return rows_sub_to_master_.at(row_in_master + size_master_ * col_in_master).at(row_in_subproblem);
}

template<class UserNodeLabel, class UserEdgeLabel>
typename BranchTight<UserNodeLabel, UserEdgeLabel>::SizeTSizeTMap_ &
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
cols_subproblem_to_master(std::size_t row_in_master, std::size_t col_in_master) {
	return cols_sub_to_master_.at(row_in_master + size_master_ * col_in_master);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
col_in_master(std::size_t row_in_master, std::size_t col_in_master, std::size_t col_in_subproblem) const {
	return cols_sub_to_master_.at(row_in_master + size_master_ * col_in_master).at(col_in_subproblem);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
get_size() const {
	return size_master_;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
get_degree() const {
	return degree_;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchTight<UserNodeLabel, UserEdgeLabel>::
SubproblemSolvers_ ::
solve(std::size_t num_threads) {
#ifdef _OPENMP
	omp_set_num_threads(num_threads - 1);
#pragma omp parallel for if(num_threads > 1)
#endif
	for (std::size_t i = 0; i < subproblem_solvers_.size(); i++) {
		subproblem_solvers_.at(i).solve();
	}
}

}

#endif /* SRC_METHODS_BRANCH_TIGHT_IPP_ */
