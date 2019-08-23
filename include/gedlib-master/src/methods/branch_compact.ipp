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
 * @file  branch_compact.ipp
 * @brief BranchCompact class definition.
 */

#ifndef SRC_METHODS_BRANCH_COMPACT_IPP_
#define SRC_METHODS_BRANCH_COMPACT_IPP_

namespace ged {

// === Definitions of destructors and constructors. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchCompact<UserNodeLabel, UserEdgeLabel>::
~BranchCompact() {}

template<class UserNodeLabel, class UserEdgeLabel>
BranchCompact<UserNodeLabel, UserEdgeLabel>::
BranchCompact(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
sort_method_{COUNTING},
branches_() {}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BranchCompact<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
	for (auto graph = this->ged_data_.begin(); graph != this->ged_data_.end(); graph++) {
		init_graph_(*graph);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchCompact<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {
	// Initialize the branches.
	if (not this->initialized_) {
		init_graph_(g);
		init_graph_(h);
	}
	std::list<Branch_> branches_g(branches_.at(g.id()));
	std::list<Branch_> branches_h(branches_.at(h.id()));
	double min_edit_cost{this->ged_data_.min_edit_cost(g, h)};
	// Delete common branches.
	auto branch_g = branches_g.begin();
	auto branch_h = branches_h.begin();
	while ((branch_g != branches_g.end()) and (branch_h != branches_h.end())) {
		int comp{branch_g->compare(*branch_h)};
		if (comp == 0) {
			branch_g = branches_g.erase(branch_g);
			branch_h = branches_h.erase(branch_h);
		}
		else if (comp == -1) {
			branch_g++;
		}
		else {
			branch_h++;
		}
	}
	// Compute lower bound from remaining branches.
	double lower_bound{0.0};
	branch_g = branches_g.begin();
	branch_h = branches_h.begin();
	while ((branch_g != branches_g.end()) and (branch_h != branches_h.end())) {
		int comp{branch_g->compare(*branch_h)};
		if (branch_g->node_label == branch_h->node_label) {
			branch_g = branches_g.erase(branch_g);
			branch_h = branches_h.erase(branch_h);
			lower_bound += 0.5 * min_edit_cost;
		}
		else if (comp == -1) {
			branch_g++;
		}
		else {
			branch_h++;
		}
	}
	lower_bound += (static_cast<double>(std::max(branches_g.size(), branches_h.size())) * min_edit_cost);
	result.set_lower_bound(lower_bound);
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchCompact<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "sort-method") {
		if (arg == "STD") {
			sort_method_ = STD;
		}
		else if (arg == "COUNTING") {
			sort_method_ = COUNTING;
		}
		else {
			throw ged::Error(std::string("Invalid argument \"") + arg  + "\" for option upper-bound. Usage: options = \"[--sort-method STD|COUNTING] [...]\"");
		}
		return true;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
BranchCompact<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	return "[--sort-method <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BranchCompact<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	sort_method_ = COUNTING;
}

// === Definition of private helper functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BranchCompact<UserNodeLabel, UserEdgeLabel>::
init_graph_(const GEDGraph & graph) {
	SortedUserEdgeLabels_ sorted_edge_labels(graph, sort_method_);
	std::list<Branch_> branches;
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		branches.emplace_back(graph.get_node_label(*node), sorted_edge_labels.get_incident_labels(*node));
	}
	branches.sort();
	branches_[graph.id()] = branches;
}

// === Definition of private class SortedUserEdgeLabels_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchCompact<UserNodeLabel, UserEdgeLabel>::
SortedUserEdgeLabels_ ::
SortedUserEdgeLabels_(const GEDGraph & g, SortMethod_ sort_method):
sorted_edge_labels_(){
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		sorted_edge_labels_[*node] = std::vector<LabelID>(0);
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
const std::vector<LabelID> &
BranchCompact<UserNodeLabel, UserEdgeLabel>::
SortedUserEdgeLabels_ ::
get_incident_labels(GEDGraph::NodeID node) const {
	return sorted_edge_labels_.at(node);
}

// === Definition of private class Branch_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BranchCompact<UserNodeLabel, UserEdgeLabel>::
Branch_ ::
Branch_(LabelID node_label, const std::vector<LabelID> & sorted_edge_labels) :
node_label{node_label},
sorted_edge_labels(sorted_edge_labels){}

template<class UserNodeLabel, class UserEdgeLabel>
BranchCompact<UserNodeLabel, UserEdgeLabel>::
Branch_ ::
Branch_(const Branch_ & branch) :
node_label{branch.node_label},
sorted_edge_labels(branch.sorted_edge_labels){}

template<class UserNodeLabel, class UserEdgeLabel>
int
BranchCompact<UserNodeLabel, UserEdgeLabel>::
Branch_ ::
compare(const Branch_ & rhs) const {
	if (node_label < rhs.node_label) {
		return -1;
	}
	if (node_label > rhs.node_label) {
		return 1;
	}
	auto label_l = sorted_edge_labels.begin();
	auto label_r = rhs.sorted_edge_labels.begin();
	while ((label_l != sorted_edge_labels.end()) and (label_r != rhs.sorted_edge_labels.end())) {
		if (*label_l == *label_r) {
			label_l++;
			label_r++;
		}
		else if (*label_l < *label_r){
			return -1;
		}
		else if (*label_l > *label_r){
			return 1;
		}
	}
	if ((label_l == sorted_edge_labels.end()) and (label_r == rhs.sorted_edge_labels.end())) {
		return 0;
	}
	if (label_l == sorted_edge_labels.end()) {
		return -1;
	}
	return 1;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchCompact<UserNodeLabel, UserEdgeLabel>::
Branch_ ::
operator<(const Branch_ & rhs) const {
	return compare(rhs) == -1;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchCompact<UserNodeLabel, UserEdgeLabel>::
Branch_ ::
operator>(const Branch_ & rhs) const {
	return compare(rhs) == 1;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BranchCompact<UserNodeLabel, UserEdgeLabel>::
Branch_ ::
operator==(const Branch_ & rhs) const {
	return compare(rhs) == 0;
}

}

#endif /* SRC_METHODS_BRANCH_COMPACT_IPP_ */
