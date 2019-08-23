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
 * @file  branch_fast.ipp
 * @brief BranchFast class definition.
 */

#ifndef SRC_METHODS_BRANCH_FAST_IPP_
#define SRC_METHODS_BRANCH_FAST_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchFast<UserNodeLabel, UserEdgeLabel>::
~BranchFast() {}

template<class UserNodeLabel, class UserEdgeLabel>
BranchFast<UserNodeLabel, UserEdgeLabel>::
BranchFast(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
sort_method_{COUNTING},
sorted_edge_labels_() {}

// === Definitions of member functions inherited from LSAPEBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BranchFast<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {
	sorted_edge_labels_[graph.id()] = SortedEdgeLabels_(graph, sort_method_);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchFast<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {
	sort_method_= COUNTING;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
BranchFast<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	return "[--sort-method <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchFast<UserNodeLabel, UserEdgeLabel>::
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
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchFast<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {

	const SortedEdgeLabels_ & sorted_edge_labels_g = sorted_edge_labels_.at(g.id());
	const SortedEdgeLabels_ & sorted_edge_labels_h = sorted_edge_labels_.at(h.id());

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < master_problem.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < master_problem.num_cols(); col_in_master++) {
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				master_problem(row_in_master, col_in_master) = compute_substitution_cost_(g, h, row_in_master, col_in_master, sorted_edge_labels_g, sorted_edge_labels_h);
			}
			else if (row_in_master < g.num_nodes()) {
				master_problem(row_in_master, h.num_nodes()) = compute_deletion_cost_(g, row_in_master);
			}
			else if (col_in_master < h.num_nodes()) {
				master_problem(g.num_nodes(), col_in_master) = compute_insertion_cost_(h, col_in_master);
			}
		}
	}

}

// === Definition of private class SortedUserEdgeLabels_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchFast<UserNodeLabel, UserEdgeLabel>::
SortedEdgeLabels_ ::
SortedEdgeLabels_(const GEDGraph & g, SortMethod_ sort_method) :
sorted_edge_labels_() {
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
BranchFast<UserNodeLabel, UserEdgeLabel>::
SortedEdgeLabels_ ::
SortedEdgeLabels_() :
sorted_edge_labels_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchFast<UserNodeLabel, UserEdgeLabel>::
SortedEdgeLabels_ ::
operator=(const SortedEdgeLabels_ & sorted_edge_labels) {
	sorted_edge_labels_ = sorted_edge_labels.sorted_edge_labels_;
}

template<class UserNodeLabel, class UserEdgeLabel>
const std::vector<LabelID> &
BranchFast<UserNodeLabel, UserEdgeLabel>::
SortedEdgeLabels_ ::
get_incident_labels(GEDGraph::NodeID node) const {
	return sorted_edge_labels_.at(node);
}

// === Definition of private helper functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
double
BranchFast<UserNodeLabel, UserEdgeLabel>::
compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k,
		const SortedEdgeLabels_ & sorted_edge_labels_g, const SortedEdgeLabels_ & sorted_edge_labels_h) const {
	// Collect node substitution cost.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k))};

	// Compute and add minimal edge insertion costs.
	if (g.degree(i) < h.degree(k)) {
		double min_edge_ins_cost{std::numeric_limits<double>::infinity()};
		for (auto label_h = sorted_edge_labels_h.get_incident_labels(k).begin(); label_h != sorted_edge_labels_h.get_incident_labels(k).end(); label_h++) {
			min_edge_ins_cost = std::min(min_edge_ins_cost, this->ged_data_.edge_cost(dummy_label(), *label_h));
		}
		cost += static_cast<double>(h.degree(k) - g.degree(i)) * min_edge_ins_cost * 0.5;
	}

	// Compute and add minimal edge deletion costs.
	if (g.degree(i) > h.degree(k)) {
		double min_edge_del_cost{std::numeric_limits<double>::infinity()};
		for (auto label_g = sorted_edge_labels_g.get_incident_labels(i).begin(); label_g != sorted_edge_labels_g.get_incident_labels(i).end(); label_g++) {
			min_edge_del_cost = std::min(min_edge_del_cost, this->ged_data_.edge_cost(*label_g, dummy_label()));
		}
		cost += static_cast<double>(g.degree(i) - h.degree(k)) * min_edge_del_cost * 0.5;
	}

	// Compute minimal edge relabelling costs.
	double min_edge_subs_cost{std::numeric_limits<double>::infinity()};
	for (auto label_g = sorted_edge_labels_g.get_incident_labels(i).begin(); label_g != sorted_edge_labels_g.get_incident_labels(i).end(); label_g++) {
		for (auto label_h = sorted_edge_labels_h.get_incident_labels(k).begin(); label_h != sorted_edge_labels_h.get_incident_labels(k).end(); label_h++) {
			if (*label_g != *label_h) {
				min_edge_subs_cost = std::min(min_edge_subs_cost, this->ged_data_.edge_cost(*label_g, *label_h));
			}
		}
	}

	// Compute multiset intersection size.
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

	// Collect edge relabelling costs.
	if (std::min(g.degree(i), h.degree(k)) - intersection_size > 0) {
		cost += static_cast<double>(std::min(g.degree(i), h.degree(k)) - intersection_size) * min_edge_subs_cost * 0.5;
	}

	// Return overall substitution cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchFast<UserNodeLabel, UserEdgeLabel>::
compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const {
	// Collect node deletion cost.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), ged::dummy_label())};

	// Collect edge deletion costs.
	auto incident_edges_i = g.incident_edges(i);
	for (auto ij = incident_edges_i.first; ij != incident_edges_i.second; ij++) {
		cost += this->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label()) * 0.5;
	}

	// Return overall deletion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
BranchFast<UserNodeLabel, UserEdgeLabel>::
compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const {
	// Collect node insertion cost.
	double cost{this->ged_data_.node_cost(ged::dummy_label(), h.get_node_label(k))};

	// Collect edge insertion costs.
	auto incident_edges_k = h.incident_edges(k);
	for (auto kl = incident_edges_k.first; kl != incident_edges_k.second; kl++) {
		cost += this->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl)) * 0.5;
	}

	// Return overall insertion cost.
	return cost;
}

}

#endif /* SRC_METHODS_BRANCH_FAST_IPP_ */
