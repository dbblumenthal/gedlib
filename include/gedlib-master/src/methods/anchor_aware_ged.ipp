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
 * @file  anchor_aware_ged.ipp
 * @brief ged::AnchorAwareGED class definition.
 */

#ifndef SRC_METHODS_ANCHOR_AWARE_GED_IPP_
#define SRC_METHODS_ANCHOR_AWARE_GED_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
~AnchorAwareGED() {}

template<class UserNodeLabel, class UserEdgeLabel>
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
AnchorAwareGED(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
lsape_model_{LSAPESolver::ECBP},
search_method_{DFS},
lower_bound_method_{BRANCH_FAST},
num_threads_{1},
time_limit_in_sec_{0.0},
map_root_to_root_{false},
sorted_edges_(),
best_feasible_(0, 0),
open_(),
omega_{0.0} {}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
	for (auto graph = this->ged_data_.begin(); graph != this->ged_data_.end(); graph++) {
		init_graph_(*graph);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {
	best_feasible_ = NodeMap(g.num_nodes(), h.num_nodes());
	open_ = std::priority_queue<TreeNode_>();
	omega_ = this->ged_data_.max_edit_cost(g, h) + 10.0;
	Timer timer(time_limit_in_sec_);

	if (not this->initialized_ and lower_bound_method_ == BRANCH_FAST) {
		init_graph_(g);
		init_graph_(h);
	}

	TreeNode_ current_node(g, h, this);
	if (current_node.is_leaf_node()) {
		current_node.extend_leaf_node(g, h);
		result.add_node_map(current_node.node_map());
		result.set_lower_bound(current_node.induced_cost());
		result.sort_node_maps_and_set_upper_bound();
		return;
	}
	generate_best_child_(g, h, current_node);

	if (not map_root_to_root_) {
		IPFP<UserNodeLabel, UserEdgeLabel> ipfp(this->ged_data_);
		ipfp.set_options("--quadratic-model QAPE --threads " + std::to_string(num_threads_));
		Result ipfp_result;
		ipfp.run_as_util(g, h, ipfp_result);
		best_feasible_ = ipfp_result.node_map(0);
	}

	while (not open_.empty() and not timer.expired()) {
		current_node = open_.top();
		open_.pop();
		if (current_node.lower_bound() >= best_feasible_.induced_cost()) {
			continue;
		}
		if (current_node.is_leaf_node()) {
			current_node.extend_leaf_node(g, h);
			if (current_node.induced_cost() < best_feasible_.induced_cost()) {
				best_feasible_ = current_node.node_map();
			}
		}
		else {
			if (current_node.has_unexplored_sibling()) {
				generate_best_sibling_(g, h, current_node);
			}
			generate_best_child_(g, h, current_node);
		}
	}

	result.add_node_map(best_feasible_);
	result.sort_node_maps_and_set_upper_bound();
	if (open_.empty()) {
		result.set_lower_bound(best_feasible_.induced_cost());
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	lsape_model_ = LSAPESolver::ECBP;
	search_method_ = DFS;
	lower_bound_method_ = BRANCH_FAST;
	num_threads_ = 1;
	time_limit_in_sec_ = 0.0;
	map_root_to_root_ = false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	return "[--lsape-model <arg>] [--search-method <arg>] [--lower-bound-method <arg>] [--threads <arg>] [--time-limit] [--map-root-to-root <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
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
		else if (arg == "FBP0") {
			lsape_model_ = LSAPESolver::FBP0;
		}
		else if (arg  == "ECBP") {
			lsape_model_ = LSAPESolver::ECBP;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option lsape-model. Usage: options = \"[--lsape-model ECBP|EBP|FLWC|FLCC|FBP|SFBP|FBP0] [...]\"");
		}
		return true;
	}
	else if (option == "search-method") {
		if (arg == "DFS") {
			search_method_ = DFS;
		}
		else if (arg  == "ASTAR") {
			search_method_ = ASTAR;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option search-method. Usage: options = \"[--search-method DFS|ASTAR] [...]\"");
		}
		return true;
	}
	else if (option == "lower-bound-method") {
		if (arg == "BRANCH") {
			lower_bound_method_ = BRANCH;
		}
		else if (arg  == "BRANCH_FAST") {
			lower_bound_method_ = BRANCH_FAST;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option lower-bound-method. Usage: options = \"[--lower-bound-method BRANCH|BRANCH_FAST] [...]\"");
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
		return true;
	}
	else if (option == "threads") {
		try {
			num_threads_ = std::stoi(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument ") + arg + " for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		if (num_threads_ <= 0) {
			throw Error(std::string("Invalid argument ") + arg + " for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		return true;
	}
	else if (option == "map-root-to-root") {
		if (arg == "TRUE") {
			map_root_to_root_ = true;
		}
		else if (arg == "FALSE") {
			map_root_to_root_ = false;
		}
		else {
			throw Error(std::string("Invalid argument ") + arg + " for option map-root-to-root. Usage: options = \"[--map-root-to-root TRUE|FALSE] [...]");
		}
		return true;
	}
	return false;
}

// ==== Definition of private struct Edge_. ====
template<class UserNodeLabel, class UserEdgeLabel>
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
Edge_ ::
Edge_(LabelID label, GEDGraph::EdgeID edge_id) :
label{label},
edge_id{edge_id}{}

template<class UserNodeLabel, class UserEdgeLabel>
bool
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
Edge_ ::
operator<(const Edge_ & rhs) const {
	return label < rhs.label;
}

// ==== Definition of private class SortedEdges_. ====
template<class UserNodeLabel, class UserEdgeLabel>
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
SortedEdges_ ::
SortedEdges_() :
sorted_edges_() {}

template<class UserNodeLabel, class UserEdgeLabel>
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
SortedEdges_ ::
SortedEdges_(const GEDGraph & g):
sorted_edges_(){
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		sorted_edges_[*node] = std::vector<Edge_>();
		for (auto edge = g.incident_edges(*node).first; edge != g.incident_edges(*node).second; edge++) {
			sorted_edges_[*node].push_back(Edge_(g.get_edge_label(*edge), *edge));
		}
		std::sort(sorted_edges_[*node].begin(), sorted_edges_[*node].end());
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
SortedEdges_ ::
operator=(const SortedEdges_ & rhs) {
	sorted_edges_ = rhs.sorted_edges_;
}

template<class UserNodeLabel, class UserEdgeLabel>
const std::vector<typename AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::Edge_> &
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
SortedEdges_::
get_incident_edges(GEDGraph::NodeID node) const {
	return sorted_edges_.at(node);
}

// ==== Definition of private class TreeNode_. ====
template<class UserNodeLabel, class UserEdgeLabel>
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
TreeNode_(const GEDGraph & g, const GEDGraph & h, const AnchorAwareGED * exact) :
exact_{exact},
node_map_(g.num_nodes(), h.num_nodes()),
is_matched_node_in_g_(g.num_nodes(), false),
is_matched_node_in_h_(h.num_nodes(), false),
is_candidate_in_h_(h.num_nodes(), true),
dummy_node_is_candidate_in_h_{true},
original_id_of_unmatched_nodes_in_h_(),
induced_cost_{0.0},
lower_bound_to_leaf_{0.0},
num_matched_nodes_in_g_{0},
num_matched_nodes_in_h_{0} {
	if (exact_->map_root_to_root_) {
		num_matched_nodes_in_g_++;
		num_matched_nodes_in_h_++;
		is_matched_node_in_g_[0] = true;
		is_matched_node_in_h_[0] = true;
		is_candidate_in_h_[0] = false;
		node_map_.add_assignment(0, 0);
		induced_cost_ = exact_->ged_data_.node_cost(g.get_node_label(0), h.get_node_label(0));
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
TreeNode_(const TreeNode_ & tree_node) :
exact_{tree_node.exact_},
node_map_(tree_node.node_map_),
is_matched_node_in_g_(tree_node.is_matched_node_in_g_),
is_matched_node_in_h_(tree_node.is_matched_node_in_h_),
is_candidate_in_h_(tree_node.is_candidate_in_h_),
dummy_node_is_candidate_in_h_(tree_node.dummy_node_is_candidate_in_h_),
original_id_of_unmatched_nodes_in_h_(tree_node.original_id_of_unmatched_nodes_in_h_),
induced_cost_{tree_node.induced_cost_},
lower_bound_to_leaf_{tree_node.lower_bound_to_leaf_},
num_matched_nodes_in_g_{tree_node.num_matched_nodes_in_g_},
num_matched_nodes_in_h_{tree_node.num_matched_nodes_in_h_} {}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
operator=(const TreeNode_ & rhs) {
	exact_ = rhs.exact_;
	node_map_ = rhs.node_map_;
	is_matched_node_in_g_ = rhs.is_matched_node_in_g_;
	is_matched_node_in_h_ = rhs.is_matched_node_in_h_;
	is_candidate_in_h_ = rhs.is_candidate_in_h_;
	dummy_node_is_candidate_in_h_ = rhs.dummy_node_is_candidate_in_h_;
	original_id_of_unmatched_nodes_in_h_ = rhs.original_id_of_unmatched_nodes_in_h_;
	induced_cost_ = rhs.induced_cost_;
	lower_bound_to_leaf_ = rhs.lower_bound_to_leaf_;
	num_matched_nodes_in_g_ = rhs.num_matched_nodes_in_g_;
	num_matched_nodes_in_h_ = rhs.num_matched_nodes_in_h_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
lower_bound() const {
	return induced_cost_ + lower_bound_to_leaf_;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
operator<(const TreeNode_ & rhs) const {
	if (exact_->search_method_ == ASTAR){
		return lower_bound() > rhs.lower_bound();
	}
	return num_matched_nodes_in_g_ < rhs.num_matched_nodes_in_g_;
}

template<class UserNodeLabel, class UserEdgeLabel>
GEDGraph::NodeID
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
next_unmatched_node_in_g() const {
	if (num_matched_nodes_in_g_ < is_matched_node_in_g_.size()) {
		return num_matched_nodes_in_g_;
	}
	return GEDGraph::dummy_node();
}

template<class UserNodeLabel, class UserEdgeLabel>
GEDGraph::NodeID
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
last_matched_node_in_g() const {
	if (num_matched_nodes_in_g_ > 0) {
		return num_matched_nodes_in_g_ - 1;
	}
	return GEDGraph::undefined_node();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
prepare_for_child_generation() {
	for (GEDGraph::NodeID node_in_h{0}; node_in_h < is_matched_node_in_h_.size(); node_in_h++) {
		is_candidate_in_h_[node_in_h] = not is_matched_node_in_h_.at(node_in_h);
	}
	dummy_node_is_candidate_in_h_ = true;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
append_extension(const GEDGraph & g, const GEDGraph & h, const NodeMap & extension) {
	std::vector<NodeMap::Assignment> assignments;
	extension.as_relation(assignments);
	for (const auto & assignment : assignments) {
		if (assignment.first != GEDGraph::dummy_node() and assignment.second != GEDGraph::dummy_node()) {
			node_map_.add_assignment(assignment.first + num_matched_nodes_in_g_, original_id_of_unmatched_nodes_in_h_.at(assignment.second));
		}
		else if (assignment.first != GEDGraph::dummy_node()) {
			node_map_.add_assignment(assignment.first + num_matched_nodes_in_g_, GEDGraph::dummy_node());
		}
		else {
			node_map_.add_assignment(GEDGraph::dummy_node(), original_id_of_unmatched_nodes_in_h_.at(assignment.second));
		}
	}
	exact_->ged_data_.compute_induced_cost(g, h, node_map_);
	induced_cost_ = node_map_.induced_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
append_next_assignment(const NodeMap & extension) {
	GEDGraph::NodeID next_node_in_g = num_matched_nodes_in_g_++;
	is_matched_node_in_g_[next_node_in_g] = true;
	GEDGraph::NodeID next_node_in_h_in_extension{extension.image(0)};
	GEDGraph::NodeID next_node_in_h{next_node_in_h_in_extension == GEDGraph::dummy_node() ? GEDGraph::dummy_node() : original_id_of_unmatched_nodes_in_h_.at(next_node_in_h_in_extension)};
	if (next_node_in_h != GEDGraph::dummy_node()) {
		num_matched_nodes_in_h_++;
		is_matched_node_in_h_[next_node_in_h] = true;
		is_candidate_in_h_[next_node_in_h] = false;
	}
	else {
		dummy_node_is_candidate_in_h_ = false;
	}
	node_map_.add_assignment(next_node_in_g, next_node_in_h);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
populate_lsape_instance(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance) {
#ifdef _OPENMP
	omp_set_num_threads(exact_->num_threads_ - 1);
#pragma omp parallel for if(exact_->num_threads_ > 1)
#endif
	for (std::size_t row = 0; row < lsape_instance.num_rows(); row++) {
		for (std::size_t col = 0; col < lsape_instance.num_cols(); col++) {
			if ((row < lsape_instance.num_rows() - 1) and (col < lsape_instance.num_cols() - 1)) {
				if ((row == 0) and (not is_candidate_in_h_.at(original_id_of_unmatched_nodes_in_h_.at(col)))) {
					lsape_instance(row, col) = exact_->omega_;
				}
				else {
					if (exact_->lower_bound_method_ == BRANCH) {
						lsape_instance(row, col) = compute_branch_substitution_cost_(g, h, row + num_matched_nodes_in_g_, original_id_of_unmatched_nodes_in_h_.at(col));
					}
					else {
						lsape_instance(row, col) = compute_branch_fast_substitution_cost_(g, h, row + num_matched_nodes_in_g_, original_id_of_unmatched_nodes_in_h_.at(col));
					}
				}
			}
			else if (row < lsape_instance.num_rows() - 1) {
				if ((row == 0) and (not dummy_node_is_candidate_in_h_)) {
					lsape_instance(row, col) = exact_->omega_;
				}
				else {
					lsape_instance(row, col) = compute_deletion_cost_(g, row + num_matched_nodes_in_g_);
				}
			}
			else if (col < lsape_instance.num_cols() - 1) {
				lsape_instance(row, col) = compute_insertion_cost_(h, original_id_of_unmatched_nodes_in_h_.at(col));
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const {
	// Collect node deletion cost.
	double cost{exact_->ged_data_.node_cost(g.get_node_label(i), ged::dummy_label())};

	// Collect edge deletion costs.
	auto incident_edges_i = g.incident_edges(i);
	for (auto ij = incident_edges_i.first; ij != incident_edges_i.second; ij++) {
		if (not is_matched_node_in_g_.at(g.head(*ij))) {
			cost += exact_->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label()) * 0.5;
		}
		else {
			cost += exact_->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label());
		}
	}

	// Return overall deletion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const {
	// Collect node insertion cost.
	double cost{exact_->ged_data_.node_cost(ged::dummy_label(), h.get_node_label(k))};

	// Collect edge insertion costs.
	auto incident_edges_k = h.incident_edges(k);
	for (auto kl = incident_edges_k.first; kl != incident_edges_k.second; kl++) {
		if (not is_matched_node_in_h_.at(h.head(*kl))) {
			cost += exact_->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl)) * 0.5;
		}
		else {
			cost += exact_->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl));
		}
	}

	// Return overall insertion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
compute_branch_fast_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const {
	// Collect node substitution costs.
	double cost{exact_->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k))};

	// Collect outer edge costs.
	std::vector<NodeMap::Assignment> assignments;
	node_map_.as_relation(assignments);
	for (const auto & assignment : assignments) {
		GEDGraph::NodeID j{assignment.first};
		GEDGraph::NodeID l{assignment.second};
		if (g.is_edge(i, j) and h.is_edge(k, l)) {
			cost += exact_->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), h.get_edge_label(h.get_edge(k, l)));
		}
		else if (g.is_edge(i, j)) {
			cost += exact_->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), dummy_label());
		}
		else if (h.is_edge(k, l)) {
			cost += exact_->ged_data_.edge_cost(dummy_label(), h.get_edge_label(h.get_edge(k, l)));
		}
	}

	// Collect unmatched edge labels.
	std::vector<LabelID> edge_labels_to_unmatched_neighbours_i;
	for (auto ij = exact_->sorted_edges_.at(g.id()).get_incident_edges(i).begin(); ij != exact_->sorted_edges_.at(g.id()).get_incident_edges(i).end(); ij++) {
		if (not is_matched_node_in_g_.at(g.head(ij->edge_id))) {
			edge_labels_to_unmatched_neighbours_i.push_back(ij->label);
		}
	}
	std::vector<LabelID> edge_labels_to_unmatched_neighbours_k;
	for (auto kl = exact_->sorted_edges_.at(h.id()).get_incident_edges(k).begin(); kl != exact_->sorted_edges_.at(h.id()).get_incident_edges(k).end(); kl++) {
		if (not is_matched_node_in_h_.at(h.head(kl->edge_id))) {
			edge_labels_to_unmatched_neighbours_k.push_back(kl->label);
		}
	}

	// Compute and add minimal edge insertion costs.
	if (edge_labels_to_unmatched_neighbours_i.size() < edge_labels_to_unmatched_neighbours_k.size()) {
		double min_ins_cost{std::numeric_limits<double>::infinity()};
		for (auto label_h = edge_labels_to_unmatched_neighbours_k.begin(); label_h != edge_labels_to_unmatched_neighbours_k.end(); label_h++) {
			min_ins_cost = std::min(min_ins_cost, exact_->ged_data_.edge_cost(dummy_label(), *label_h));
		}
		cost += static_cast<double>(edge_labels_to_unmatched_neighbours_k.size() - edge_labels_to_unmatched_neighbours_i.size()) * min_ins_cost * 0.5;
	}

	// Compute and add minimal edge deletion costs.
	if (edge_labels_to_unmatched_neighbours_i.size() > edge_labels_to_unmatched_neighbours_k.size()) {
		double min_del_cost{std::numeric_limits<double>::infinity()};
		for (auto label_g = edge_labels_to_unmatched_neighbours_i.begin(); label_g != edge_labels_to_unmatched_neighbours_i.end(); label_g++) {
			min_del_cost = std::min(min_del_cost, exact_->ged_data_.edge_cost(*label_g, dummy_label()));
		}
		cost += static_cast<double>(edge_labels_to_unmatched_neighbours_i.size() - edge_labels_to_unmatched_neighbours_k.size()) * min_del_cost * 0.5;
	}

	// Compute minimal edge relabelling costs.
	double min_rel_cost{std::numeric_limits<double>::infinity()};
	for (auto label_g = edge_labels_to_unmatched_neighbours_i.begin(); label_g != edge_labels_to_unmatched_neighbours_i.end(); label_g++) {
		for (auto label_h = edge_labels_to_unmatched_neighbours_k.begin(); label_h != edge_labels_to_unmatched_neighbours_k.end(); label_h++) {
			if (*label_g != *label_h) {
				min_rel_cost = std::min(min_rel_cost, exact_->ged_data_.edge_cost(*label_g, *label_h));
			}
		}
	}

	// Compute multiset intersection size.
	std::size_t intersection_size{0};
	auto label_g = edge_labels_to_unmatched_neighbours_i.begin();
	auto label_h = edge_labels_to_unmatched_neighbours_k.begin();
	while ((label_g != edge_labels_to_unmatched_neighbours_i.end()) and (label_h != edge_labels_to_unmatched_neighbours_k.end())) {
		if (*label_g == *label_h) {
			intersection_size++;
			label_g++;
			label_h++;
		}
		else if (*label_g < *label_h) {
			label_g++;
		}
		else {
			label_h++;
		}
	}

	// Collect edge relabelling costs.
	std::size_t gamma(std::min(edge_labels_to_unmatched_neighbours_i.size(), edge_labels_to_unmatched_neighbours_k.size()) - intersection_size);
	if (gamma > 0) {
		cost += static_cast<double>(gamma) * min_rel_cost * 0.5;
	}

	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
compute_branch_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const {
	// Collect node substitution costs.
	double cost{exact_->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k))};

	// Collect outer edge costs.
	std::vector<NodeMap::Assignment> assignments;
	node_map_.as_relation(assignments);
	for (const auto & assignment : assignments) {
		GEDGraph::NodeID j{assignment.first};
		GEDGraph::NodeID l{assignment.second};
		if (g.is_edge(i, j) and h.is_edge(k, l)) {
			cost += exact_->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), h.get_edge_label(h.get_edge(k, l)));
		}
		else if (g.is_edge(i, j)) {
			cost += exact_->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), dummy_label());
		}
		else if (h.is_edge(k, l)) {
			cost += exact_->ged_data_.edge_cost(dummy_label(), h.get_edge_label(h.get_edge(k, l)));
		}
	}

	// Initialize subproblem.
	std::vector<LabelID> edge_labels_to_unmatched_neighbours_i;
	for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++) {
		if (not is_matched_node_in_g_.at(g.head(*ij))) {
			edge_labels_to_unmatched_neighbours_i.push_back(g.get_edge_label(*ij));
		}
	}
	std::vector<LabelID> edge_labels_to_unmatched_neighbours_k;
	for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++) {
		if (not is_matched_node_in_h_.at(h.head(*kl))) {
			edge_labels_to_unmatched_neighbours_k.push_back(h.get_edge_label(*kl));
		}
	}
	DMatrix subproblem(edge_labels_to_unmatched_neighbours_i.size() + 1, edge_labels_to_unmatched_neighbours_k.size() + 1);

	// Collect edge deletion costs.
	std::size_t row{0};
	for (auto label_ij = edge_labels_to_unmatched_neighbours_i.begin(); label_ij != edge_labels_to_unmatched_neighbours_i.end(); label_ij++, row++) {
		subproblem(row, edge_labels_to_unmatched_neighbours_k.size()) = exact_->ged_data_.edge_cost(*label_ij, ged::dummy_label()) * 0.5;
	}

	// Collect edge insertion costs.
	std::size_t col{0};
	for (auto label_kl = edge_labels_to_unmatched_neighbours_k.begin(); label_kl != edge_labels_to_unmatched_neighbours_k.end(); label_kl++, col++) {
		subproblem(edge_labels_to_unmatched_neighbours_i.size(), col) = exact_->ged_data_.edge_cost(ged::dummy_label(), *label_kl) * 0.5;
	}

	// Collect edge relabelling costs.
	row = 0;
	for (auto label_ij = edge_labels_to_unmatched_neighbours_i.begin(); label_ij != edge_labels_to_unmatched_neighbours_i.end(); label_ij++, row++) {
		col = 0;
		for (auto label_kl = edge_labels_to_unmatched_neighbours_k.begin(); label_kl != edge_labels_to_unmatched_neighbours_k.end(); label_kl++, col++) {
			subproblem(row, col) = exact_->ged_data_.edge_cost(*label_ij, *label_kl) * 0.5;
		}
	}

	// Solve subproblem.
	LSAPESolver subproblem_solver(&subproblem);
	subproblem_solver.set_model(exact_->lsape_model_);
	subproblem_solver.solve();

	// Update and return overall substitution cost.
	cost += subproblem_solver.minimal_cost();
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
set_lower_bound_to_leaf(double lower_bound_to_leaf) {
	lower_bound_to_leaf_ = lower_bound_to_leaf;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
update_induced_cost(const GEDGraph & g, const GEDGraph & h) {
	GEDGraph::NodeID i{last_matched_node_in_g()};
	GEDGraph::NodeID k{node_map_.image(i)};
	if (k != GEDGraph::dummy_node()) {
		induced_cost_ += exact_->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k));
	}
	else {
		induced_cost_ += exact_->ged_data_.node_cost(g.get_node_label(i), dummy_label());
	}
	std::vector<NodeMap::Assignment> assignments;
	node_map_.as_relation(assignments);
	for (const auto & assignment : assignments) {
		GEDGraph::NodeID j{assignment.first};
		GEDGraph::NodeID l{assignment.second};
		if (g.is_edge(i, j) and h.is_edge(k, l)) {
			induced_cost_ += exact_->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), h.get_edge_label(h.get_edge(k, l)));
		}
		else if (g.is_edge(i, j)) {
			induced_cost_ += exact_->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), dummy_label());
		}
		else if (h.is_edge(k, l)) {
			induced_cost_ += exact_->ged_data_.edge_cost(dummy_label(), h.get_edge_label(h.get_edge(k, l)));
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
induced_cost() const {
	return induced_cost_;
}

template<class UserNodeLabel, class UserEdgeLabel>
const NodeMap &
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
node_map() const {
	return node_map_;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
is_leaf_node() const {
	return ((num_matched_nodes_in_g_ == is_matched_node_in_g_.size()) or (num_matched_nodes_in_h_ == is_matched_node_in_h_.size()));
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
extend_leaf_node(const GEDGraph & g, const GEDGraph & h) {
	for (GEDGraph::NodeID i{0}; i < is_matched_node_in_g_.size(); i++) {
		if (not is_matched_node_in_g_.at(i)) {
			node_map_.add_assignment(i, GEDGraph::dummy_node());
			is_matched_node_in_g_[i] = true;
			num_matched_nodes_in_g_++;
		}
	}
	for (GEDGraph::NodeID k{0}; k < is_matched_node_in_h_.size(); k++) {
		if (not is_matched_node_in_h_.at(k)) {
			node_map_.add_assignment(GEDGraph::dummy_node(), k);
			is_matched_node_in_h_[k] = true;
			num_matched_nodes_in_h_++;
		}
	}
	exact_->ged_data_.compute_induced_cost(g, h, node_map_);
	induced_cost_ = node_map_.induced_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
prepare_for_sibling_generation() {
	GEDGraph::NodeID next_node_g{--num_matched_nodes_in_g_};
	is_matched_node_in_g_[next_node_g] = false;
	GEDGraph::NodeID next_node_h{node_map_.image(next_node_g)};
	node_map_.erase_image(next_node_g);
	if (next_node_h != GEDGraph::dummy_node()) {
		is_matched_node_in_h_[next_node_h] = false;
		num_matched_nodes_in_h_--;
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
has_unexplored_sibling() {
	if (num_matched_nodes_in_g_ == 0) {
		return false;
	}
	if (dummy_node_is_candidate_in_h_) {
		return true;
	}
	for (auto is_candidate : is_candidate_in_h_) {
		if (is_candidate) {
			return true;
		}
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
num_unmatched_nodes_in_g() const {
	return is_matched_node_in_g_.size() - num_matched_nodes_in_g_;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
num_unmatched_nodes_in_h() const {
	return is_matched_node_in_h_.size() - num_matched_nodes_in_h_;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
TreeNode_::
update_original_id_of_unmatched_nodes_in_h() {
	original_id_of_unmatched_nodes_in_h_.clear();
	GEDGraph::NodeID k{0};
	for (auto is_matched_node : is_matched_node_in_h_) {
		if (not is_matched_node) {
			original_id_of_unmatched_nodes_in_h_.push_back(k);
		}
	}
}

// ==== Definitions of private helper member functions. ====
template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
init_graph_(const GEDGraph & graph) {
	sorted_edges_[graph.id()] = SortedEdges_(graph);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
generate_next_tree_node_(const GEDGraph & g, const GEDGraph & h, TreeNode_ & next_tree_node, bool update_induced_cost, bool update_upper_bound) {

	next_tree_node.update_original_id_of_unmatched_nodes_in_h();
	// construct LSAPE instance
	DMatrix lsape_instance(next_tree_node.num_unmatched_nodes_in_g() + 1, next_tree_node.num_unmatched_nodes_in_h() + 1, 0.0);
	next_tree_node.populate_lsape_instance(g, h, lsape_instance);

	// solve LSAPE instance and update lower bound to leaf
	LSAPESolver lsape_solver(&lsape_instance);
	lsape_solver.set_model(lsape_model_);
	lsape_solver.solve();
	next_tree_node.set_lower_bound_to_leaf(lsape_solver.minimal_cost());

	// update node map
	NodeMap extension(next_tree_node.num_unmatched_nodes_in_g(), next_tree_node.num_unmatched_nodes_in_h());
	util::construct_node_map_from_solver(lsape_solver, extension);
	next_tree_node.append_next_assignment(extension);

	// update induced cost
	if (update_induced_cost) {
		next_tree_node.update_induced_cost(g, h);
	}
	open_.push(next_tree_node);

	if (update_upper_bound) {
		extension.erase_image(0);
		next_tree_node.append_extension(g, h, extension);
		if (next_tree_node.induced_cost() < best_feasible_.induced_cost()) {
			best_feasible_ = next_tree_node.node_map();
		}
	}

}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
generate_best_child_(const GEDGraph & g, const GEDGraph & h, const TreeNode_ & current_node) {
	TreeNode_ child_node(current_node);
	child_node.prepare_for_child_generation();
	generate_next_tree_node_(g, h, child_node, true, false);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::
generate_best_sibling_(const GEDGraph & g, const GEDGraph & h, const TreeNode_ & current_node) {
	TreeNode_ sibling_node(current_node);
	sibling_node.prepare_for_sibling_generation();
	generate_next_tree_node_(g, h, sibling_node, false, false);
}

}

#endif /* SRC_METHODS_ANCHOR_AWARE_GED_IPP_ */
