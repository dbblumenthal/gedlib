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
 * @file  branch_uniform.ipp
 * @brief BranchUniform class definition.
 */

#ifndef SRC_METHODS_BRANCH_UNIFORM_IPP_
#define SRC_METHODS_BRANCH_UNIFORM_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchUniform<UserNodeLabel, UserEdgeLabel>::
~BranchUniform() {}

template<class UserNodeLabel, class UserEdgeLabel>
BranchUniform<UserNodeLabel, UserEdgeLabel>::
BranchUniform(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
sort_method_{COUNTING},
wildcard_option_{false},
sorted_edge_labels_() {}

// === Definitions of member functions inherited from LSAPEBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {
	sorted_edge_labels_[graph.id()] = SortedEdgeLabels_(graph, sort_method_);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {

	const SortedEdgeLabels_ & sorted_edge_labels_g = sorted_edge_labels_.at(g.id());
	const SortedEdgeLabels_ & sorted_edge_labels_h = sorted_edge_labels_.at(h.id());
	double min_edge_subs_cost{this->ged_data_.min_edge_subs_cost(g, h)};
	double min_edge_del_cost{this->ged_data_.min_edge_del_cost(g)};
	double min_edge_ins_cost{this->ged_data_.min_edge_ins_cost(h)};

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < master_problem.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < master_problem.num_cols(); col_in_master++) {
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				master_problem(row_in_master, col_in_master) = compute_substitution_cost_(g, h, row_in_master, col_in_master, sorted_edge_labels_g, sorted_edge_labels_h, min_edge_subs_cost, min_edge_del_cost, min_edge_ins_cost);
			}
			else if (row_in_master < g.num_nodes()) {
				master_problem(row_in_master, h.num_nodes()) = compute_deletion_cost_(g, row_in_master, min_edge_del_cost);
			}
			else if (col_in_master < h.num_nodes()) {
				master_problem(g.num_nodes(), col_in_master) = compute_insertion_cost_(h, col_in_master, min_edge_ins_cost);
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {
	sort_method_ = COUNTING;
	wildcard_option_ = false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	return "[--sort-method <arg>] [--wildcards <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchUniform<UserNodeLabel, UserEdgeLabel>::
lsape_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "sort-method") {
		if (arg == "STD") {
			sort_method_ = STD;
		}
		else if (arg == "COUNTING") {
			sort_method_ = COUNTING;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option upper-bound. Usage: options = \"[--sort-method STD|COUNTING] [...]\"");
		}
		return true;
	}
	else if (option == "wildcards") {
		if (arg == "NO") {
			wildcard_option_ = false;
		}
		else if (arg == "YES") {
			wildcard_option_ = true;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option wildcards. Usage: options = \"[--wildcards YES|NO] [...]\"");
		}
		return true;
	}
	return false;
}

// === Definition of private class SortedUserEdgeLabels_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchUniform<UserNodeLabel, UserEdgeLabel>::
SortedEdgeLabels_ ::
SortedEdgeLabels_(const GEDGraph & g, SortMethod_ sort_method):
sorted_edge_labels_(){
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		sorted_edge_labels_[*node] = std::vector<LabelID>();
		for (auto edge = g.incident_edges(*node).first; edge != g.incident_edges(*node).second; edge++) {
			sorted_edge_labels_[*node].push_back(g.get_edge_label(*edge));
		}
		switch (sort_method) {
		case STD:
			std::sort(sorted_edge_labels_[*node].begin(), sorted_edge_labels_[*node].end());
			break;
		default:
			util::counting_sort(sorted_edge_labels_[*node].begin(), sorted_edge_labels_[*node].end());
			break;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
BranchUniform<UserNodeLabel, UserEdgeLabel>::
SortedEdgeLabels_ ::
SortedEdgeLabels_():
sorted_edge_labels_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchUniform<UserNodeLabel, UserEdgeLabel>::
SortedEdgeLabels_ ::
operator=(const SortedEdgeLabels_ & sorted_edge_labels) {
	sorted_edge_labels_ = sorted_edge_labels.sorted_edge_labels_;
}

template<class UserNodeLabel, class UserEdgeLabel>
const std::vector<LabelID> &
BranchUniform<UserNodeLabel, UserEdgeLabel>::
SortedEdgeLabels_ ::
get_incident_labels(GEDGraph::NodeID node) const {
	return sorted_edge_labels_.at(node);
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
double
BranchUniform<UserNodeLabel, UserEdgeLabel>::
compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k,
		const SortedEdgeLabels_ & sorted_edge_labels_g, const SortedEdgeLabels_ & sorted_edge_labels_h,
		double min_edge_subs_cost, double min_edge_del_cost, double min_edge_ins_cost) const {

	double cost{0.0};

	// Collect node substitution cost.
	if ((not wildcard_option_) or (h.get_node_label(k) != dummy_label())) {
		cost += this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k));
	}

	// Determine the number of wildcard edges.
	std::size_t num_incident_wildard_edges{0};
	if (wildcard_option_) {
		for (auto label_h : sorted_edge_labels_h.get_incident_labels(k)) {
			if (label_h == dummy_label()) {
				num_incident_wildard_edges++;
			}
		}
	}

	// Compute the size of the multiset intersection.
	std::size_t intersection_size{0};
	auto label_g = sorted_edge_labels_g.get_incident_labels(i).begin();
	auto label_h = sorted_edge_labels_h.get_incident_labels(k).begin();
	while ((label_g != sorted_edge_labels_g.get_incident_labels(i).end()) and (label_h != sorted_edge_labels_h.get_incident_labels(k).end())) {
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

	// Add edge deletion costs.
	if (g.degree(i) > h.degree(k)) {
		cost += static_cast<double>(g.degree(i) - h.degree(k)) * min_edge_del_cost * 0.5;
	}

	// Add edge insertion costs.
	if (g.degree(i) < h.degree(k)) {
		std::size_t num_inserted_edges{h.degree(k) - g.degree(i)};
		// Use as many wildcard edges as possible for insertion, if insertion is at least as expensive as substitution.
		if (wildcard_option_) {
			if (min_edge_ins_cost >= min_edge_subs_cost) {
				if (num_inserted_edges >= num_incident_wildard_edges) {
					num_inserted_edges -= num_incident_wildard_edges;
					num_incident_wildard_edges = 0;
				}
				else {
					num_incident_wildard_edges -= num_inserted_edges;
					num_inserted_edges = 0;
				}
			}
			else {
				std::size_t num_substituted_edges{g.degree(i) - intersection_size};
				// Use all wildcard edges for insertion that are not needed for substitution.
				if (num_substituted_edges < num_incident_wildard_edges) {
					num_inserted_edges -= (num_incident_wildard_edges - num_substituted_edges);
					num_incident_wildard_edges = num_substituted_edges;
				}
			}
		}
		cost += static_cast<double>(num_inserted_edges) * min_edge_ins_cost * 0.5;
	}


	// Add edge substitution costs.
	cost += static_cast<double>(std::min(g.degree(i), h.degree(k)) - intersection_size - num_incident_wildard_edges) * 0.5 * min_edge_subs_cost;

	// Return the overall substitution cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchUniform<UserNodeLabel, UserEdgeLabel>::
compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i, double min_edge_del_cost) const {
	// Collect node deletion cost.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), dummy_label())};

	// Collect edge deletion cost.
	cost += static_cast<double>(g.degree(i)) * 0.5 * min_edge_del_cost;

	// Return overall deletion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchUniform<UserNodeLabel, UserEdgeLabel>::
compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k, double min_edge_ins_cost) const {
	// Collect node insertion cost.
	double cost{this->ged_data_.node_cost(dummy_label(), h.get_node_label(k))};

	// Collect edge insertion cost.
	std::size_t num_incident_wildard_edges{0};
	if (wildcard_option_) {
		for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++) {
			if (h.get_edge_label(*kl) == dummy_label()) {
				num_incident_wildard_edges++;
			}
		}
	}
	cost += static_cast<double>(h.degree(k) - num_incident_wildard_edges) * 0.5 * min_edge_ins_cost;

	// Return overall insertion cost.
	return cost;
}

}

#endif /* SRC_METHODS_BRANCH_UNIFORM_IPP_ */
