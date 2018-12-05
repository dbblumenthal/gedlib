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
 * @file star.ipp
 * @brief ged::Star class definition.
 */

#ifndef SRC_METHODS_STAR_IPP_
#define SRC_METHODS_STAR_IPP_


namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
Star<UserNodeLabel, UserEdgeLabel>::
~Star() {}

template<class UserNodeLabel, class UserEdgeLabel>
Star<UserNodeLabel, UserEdgeLabel>::
Star(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
sort_method_{COUNTING},
sorted_node_labels_() {}

// === Definitions of member functions inherited from LSAPEBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Star<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {
	sorted_node_labels_[graph.id()] = SortedNodeLabels_(graph, sort_method_);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Star<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {
	sort_method_= COUNTING;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Star<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	return "[--sort-method <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Star<UserNodeLabel, UserEdgeLabel>::
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
Star<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {

	const SortedNodeLabels_ & sorted_node_labels_g = sorted_node_labels_.at(g.id());
	const SortedNodeLabels_ & sorted_node_labels_h = sorted_node_labels_.at(h.id());

	double min_edit_cost{this->ged_data_.min_edit_cost(g, h)};

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < master_problem.num_rows(); row_in_master++) {
		for (std::size_t col_in_master = 0; col_in_master < master_problem.num_cols(); col_in_master++) {
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				master_problem(row_in_master, col_in_master) = compute_substitution_cost_(g, h, row_in_master, col_in_master, sorted_node_labels_g, sorted_node_labels_h) * min_edit_cost;
			}
			else if (row_in_master < g.num_nodes()) {
				master_problem(row_in_master, h.num_nodes()) = compute_deletion_cost_(g, row_in_master) * min_edit_cost;
			}
			else if (col_in_master < h.num_nodes()) {
				master_problem(g.num_nodes(), col_in_master) = compute_insertion_cost_(h, col_in_master) * min_edit_cost;
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Star<UserNodeLabel, UserEdgeLabel>::
lsape_lower_bound_scaling_factor_(const GEDGraph & g, const GEDGraph & h) {
	std::size_t inverse_scaling_factor{std::max(g.maxdeg(), h.maxdeg()) + 1};
	if (inverse_scaling_factor < 4) {
		inverse_scaling_factor = 4;
	}
	return 1 / static_cast<double>(inverse_scaling_factor);
}

// === Definition of private class SortedUserEdgeLabels_. ===
template<class UserNodeLabel, class UserEdgeLabel>
Star<UserNodeLabel, UserEdgeLabel>::
SortedNodeLabels_ ::
SortedNodeLabels_(const GEDGraph & g, SortMethod_ sort_method) :
sorted_node_labels_() {
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		sorted_node_labels_[*node] = std::vector<LabelID>();
		for (auto edge = g.incident_edges(*node).first; edge != g.incident_edges(*node).second; edge++) {
			sorted_node_labels_[*node].push_back(g.get_node_label(g.head(*edge)));
		}
		switch (sort_method) {
		case STD:
			std::sort(sorted_node_labels_[*node].begin(), sorted_node_labels_[*node].end());
			break;
		default:
			util::counting_sort(sorted_node_labels_[*node].begin(), sorted_node_labels_[*node].end());
			break;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
Star<UserNodeLabel, UserEdgeLabel>::
SortedNodeLabels_ ::
SortedNodeLabels_() :
sorted_node_labels_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
Star<UserNodeLabel, UserEdgeLabel>::
SortedNodeLabels_ ::
operator=(const SortedNodeLabels_ & sorted_edge_labels) {
	sorted_node_labels_ = sorted_edge_labels.sorted_node_labels_;
}

template<class UserNodeLabel, class UserEdgeLabel>
const std::vector<LabelID> &
Star<UserNodeLabel, UserEdgeLabel>::
SortedNodeLabels_ ::
get_incident_labels(GEDGraph::NodeID node) const {
	return sorted_node_labels_.at(node);
}

// === Definition of private helper functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
double
Star<UserNodeLabel, UserEdgeLabel>::
compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k,
		const SortedNodeLabels_ & sorted_node_labels_g, const SortedNodeLabels_ & sorted_node_labels_h) const {

	// Compute multiset intersection size.
	std::size_t intersection_size{0};
	auto label_g = sorted_node_labels_g.get_incident_labels(i).begin();
	auto label_h = sorted_node_labels_h.get_incident_labels(k).begin();
	while ((label_g != sorted_node_labels_g.get_incident_labels(i).end()) and (label_h != sorted_node_labels_h.get_incident_labels(k).end())) {
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

	// Collect node cost.
	double cost{g.get_node_label(i) != h.get_node_label(k) ? 1.0 : 0.0};

	// Collect edge and neighbor costs.
	cost += static_cast<double>(2 * std::max(g.degree(i), h.degree(k)) - std::min(g.degree(i), h.degree(k)) - intersection_size);

	// Return overall substitution cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Star<UserNodeLabel, UserEdgeLabel>::
compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const {
	return static_cast<double>(1 + 2 * g.degree(i));
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Star<UserNodeLabel, UserEdgeLabel>::
compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const {
	return static_cast<double>(1 + 2 * h.degree(k));
}

}



#endif /* SRC_METHODS_STAR_IPP_ */
