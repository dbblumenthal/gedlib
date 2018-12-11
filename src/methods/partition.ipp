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
 * @file  partition.ipp
 * @brief Partition class definition.
 */

#ifndef SRC_METHODS_PARTITION_IPP_
#define SRC_METHODS_PARTITION_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
Partition<UserNodeLabel, UserEdgeLabel>::
~Partition() {}

template<class UserNodeLabel, class UserEdgeLabel>
Partition<UserNodeLabel, UserEdgeLabel>::
Partition(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
substruct_maps_(),
unmatched_substructs_() {}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
	for (auto graph = this->ged_data_.begin(); graph != this->ged_data_.end(); graph++) {
		init_graph_(*graph);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {

	unmatched_substructs_.clear();
	init_graphs_(g, h);

	SubstructMap_ & is_substruct_in_g = substruct_maps_.at(g.id());
	unmatched_substructs_.clear();

	std::map<GEDGraph::NodeID, bool> is_deleted_node;
	for (auto node_1 = h.nodes().first; node_1 != h.nodes().second; node_1++) {
		is_deleted_node[*node_1] = false;
	}
	std::map<GEDGraph::EdgeID, bool> is_deleted_edge;
	for (auto edge_1 = h.edges().first; edge_1 != h.edges().second; edge_1++) {
		is_deleted_edge[*edge_1] = false;
	}

	check_node_subtructs_(h, is_substruct_in_g, is_deleted_node);
	check_edge_subtructs_(h, is_substruct_in_g, is_deleted_edge);
	check_node_edge_subtructs_(h, is_substruct_in_g, is_deleted_node, is_deleted_edge);
	check_node_edge_node_subtructs_(h, is_substruct_in_g, is_deleted_node, is_deleted_edge);
	check_node_edge_edge_subtructs_(h, is_substruct_in_g, is_deleted_node, is_deleted_edge);

	result.set_lower_bound(static_cast<double>(unmatched_substructs_.size()) * this->ged_data_.min_edit_cost(g, h));
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
init_graph_(const GEDGraph & graph) {
	substruct_maps_[graph.id()] = SubstructMap_();
	SubstructMap_ & is_substruct = substruct_maps_.at(graph.id());
	for (LabelID edge_label_1{1}; edge_label_1 < this->ged_data_.num_edge_labels(); edge_label_1++) {
		Substruct_ edge_substruct(EDGE, edge_label_1);
		is_substruct[edge_substruct] = contains_substruct_(graph, edge_substruct);
	}
	for (LabelID node_label_1{1}; node_label_1 < this->ged_data_.num_node_labels(); node_label_1++) {
		Substruct_ node_substruct(NODE, node_label_1);
		bool contains_node_substruct{contains_substruct_(graph, node_substruct)};
		is_substruct[node_substruct] = contains_node_substruct;
		for (LabelID edge_label_1{1}; edge_label_1 < this->ged_data_.num_edge_labels(); edge_label_1++) {
			Substruct_ node_edge_substruct(NODE_EDGE, node_label_1, edge_label_1);
			bool contains_node_edge_substruct{false};
			if (contains_node_substruct) {
				contains_node_edge_substruct = contains_substruct_(graph, node_edge_substruct);
				is_substruct[node_edge_substruct] = contains_node_edge_substruct;
			}
			else {
				is_substruct[node_edge_substruct] = false;
			}
			for (LabelID node_label_2{1}; node_label_2 < this->ged_data_.num_node_labels(); node_label_2++) {
				Substruct_ node_edge_node_substruct(NODE_EDGE_NODE, node_label_1, edge_label_1, node_label_2);
				if (contains_node_edge_substruct) {
					is_substruct[node_edge_node_substruct] = contains_substruct_(graph, node_edge_node_substruct);
				}
				else {
					is_substruct[node_edge_node_substruct] = false;
				}
			}
			for (LabelID edge_label_2{1}; edge_label_2 < this->ged_data_.num_edge_labels(); edge_label_2++) {
				Substruct_ node_edge_edge_substruct(NODE_EDGE_EDGE, node_label_1, edge_label_1, edge_label_2);
				if (contains_node_edge_substruct) {
					is_substruct[node_edge_edge_substruct] = contains_substruct_(graph, node_edge_edge_substruct);
				}
				else {
					is_substruct[node_edge_edge_substruct] = false;
				}
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
init_graphs_(const GEDGraph & g, const GEDGraph & h) {

	std::vector<LabelID> node_labels_h;
	for (auto node = h.nodes().first; node != h.nodes().second; node++) {
		node_labels_h.push_back(h.get_node_label(*node));
	}

	std::vector<LabelID> edge_labels_h;
	for (auto edge = h.edges().first; edge != h.edges().second; edge++) {
		edge_labels_h.push_back(h.get_edge_label(*edge));
	}

	substruct_maps_[g.id()] = SubstructMap_();
	SubstructMap_ & is_substruct = substruct_maps_.at(g.id());
	for (LabelID edge_label_1 : edge_labels_h) {
		Substruct_ edge_substruct(EDGE, edge_label_1);
		is_substruct[edge_substruct] = contains_substruct_(g, edge_substruct);
	}
	for (LabelID node_label_1 : node_labels_h) {
		Substruct_ node_substruct(NODE, node_label_1);
		bool contains_node_substruct{contains_substruct_(g, node_substruct)};
		is_substruct[node_substruct] = contains_node_substruct;
		for (LabelID edge_label_1 : edge_labels_h) {
			Substruct_ node_edge_substruct(NODE_EDGE, node_label_1, edge_label_1);
			bool contains_node_edge_substruct{false};
			if (contains_node_substruct) {
				contains_node_edge_substruct = contains_substruct_(g, node_edge_substruct);
				is_substruct[node_edge_substruct] = contains_node_edge_substruct;
			}
			else {
				is_substruct[node_edge_substruct] = false;
			}
			for (LabelID node_label_2 : node_labels_h) {
				Substruct_ node_edge_node_substruct(NODE_EDGE_NODE, node_label_1, edge_label_1, node_label_2);
				if (contains_node_edge_substruct) {
					is_substruct[node_edge_node_substruct] = contains_substruct_(g, node_edge_node_substruct);
				}
				else {
					is_substruct[node_edge_node_substruct] = false;
				}
			}
			for (LabelID edge_label_2 : edge_labels_h) {
				Substruct_ node_edge_edge_substruct(NODE_EDGE_EDGE, node_label_1, edge_label_1, edge_label_2);
				if (contains_node_edge_substruct) {
					is_substruct[node_edge_edge_substruct] = contains_substruct_(g, node_edge_edge_substruct);
				}
				else {
					is_substruct[node_edge_edge_substruct] = false;
				}
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Partition<UserNodeLabel, UserEdgeLabel>::
contains_substruct_(const GEDGraph & g, const Substruct_ & substruct) const {
	switch (substruct.type) {
	case NODE:
		return contains_node_substruct_(g, substruct);
	case EDGE:
		return contains_edge_substruct_(g, substruct);
	case NODE_EDGE:
		return contains_node_edge_substruct_(g, substruct);
	case NODE_EDGE_NODE:
		return contains_node_edge_node_substruct_(g, substruct);
	case NODE_EDGE_EDGE:
		return contains_node_edge_edge_substruct_(g, substruct);
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
check_node_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::NodeID, bool> & is_deleted_node) {
	for (auto node_1 = h.nodes().first; node_1 != h.nodes().second; node_1++) {
		Substruct_ node_substruct(NODE, h.get_node_label(*node_1), dummy_label(), dummy_label(), *node_1);
		if (not is_substruct_in_g.at(node_substruct)) {
			unmatched_substructs_.push_back(node_substruct);
			is_deleted_node.at(*node_1) = true;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Partition<UserNodeLabel, UserEdgeLabel>::
contains_node_substruct_(const GEDGraph & g, const Substruct_ & substruct) const {
	for (auto node_1 = g.nodes().first; node_1 != g.nodes().second; node_1++) {
		if (g.get_node_label(*node_1) == substruct.node_label_1) {
			return true;
		}
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
check_edge_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::EdgeID, bool> & is_deleted_edge) {
	for (auto edge_1 = h.edges().first; edge_1 != h.edges().second; edge_1++) {
		Substruct_ edge_substruct(EDGE, h.get_edge_label(*edge_1), dummy_label(), dummy_label(), GEDGraph::dummy_node(), *edge_1);
		if (not is_substruct_in_g.at(edge_substruct)) {
			unmatched_substructs_.push_back(edge_substruct);
			is_deleted_edge.at(*edge_1) = true;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Partition<UserNodeLabel, UserEdgeLabel>::
contains_edge_substruct_(const GEDGraph & g, const Substruct_ & substruct) const {
	for (auto edge_1 = g.edges().first; edge_1 != g.edges().second; edge_1++) {
		if (g.get_edge_label(*edge_1) == substruct.edge_label_1) {
			return true;
		}
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
check_node_edge_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::NodeID, bool> & is_deleted_node, std::map<GEDGraph::EdgeID, bool> & is_deleted_edge) {
	for (auto node_1 = h.nodes().first; node_1 != h.nodes().second; node_1++) {
		if (is_deleted_node.at(*node_1)) {
			continue;
		}
		for (auto edge_1 = h.incident_edges(*node_1).first; edge_1 != h.incident_edges(*node_1).second; edge_1++) {
			if (is_deleted_edge.at(*edge_1)) {
				continue;
			}
			Substruct_ node_edge_substruct(NODE_EDGE, h.get_node_label(*node_1), h.get_edge_label(*edge_1), dummy_label(), *node_1, *edge_1);
			if (not is_substruct_in_g.at(node_edge_substruct)) {
				unmatched_substructs_.push_back(node_edge_substruct);
				is_deleted_node.at(*node_1) = true;
				is_deleted_edge.at(*edge_1) = true;
				break;
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Partition<UserNodeLabel, UserEdgeLabel>::
contains_node_edge_substruct_(const GEDGraph & g, const Substruct_ & substruct) const {
	std::vector<GEDGraph::NodeID> candidate_nodes;
	for (auto node_1 = g.nodes().first; node_1 != g.nodes().second; node_1++) {
		if (g.get_node_label(*node_1) == substruct.node_label_1) {
			candidate_nodes.push_back(*node_1);
		}
	}
	if (candidate_nodes.empty()) {
		return false;
	}
	for (auto node_1 = candidate_nodes.begin(); node_1 != candidate_nodes.end(); node_1++) {
		for (auto edge_1 = g.incident_edges(*node_1).first; edge_1 != g.incident_edges(*node_1).second; edge_1++) {
			if (g.get_edge_label(*edge_1) == substruct.edge_label_1) {
				return true;
			}
		}
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
check_node_edge_node_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::NodeID, bool> & is_deleted_node, std::map<GEDGraph::EdgeID, bool> & is_deleted_edge) {
	for (auto node_1 = h.nodes().first; node_1 != h.nodes().second; node_1++) {
		if (is_deleted_node.at(*node_1)) {
			continue;
		}
		for (auto edge_1 = h.incident_edges(*node_1).first; edge_1 != h.incident_edges(*node_1).second; edge_1++) {
			GEDGraph::NodeID node_2{h.head(*edge_1)};
			if (is_deleted_edge.at(*edge_1) or is_deleted_node.at(node_2)) {
				continue;
			}
			Substruct_ node_edge_node_substruct(NODE_EDGE_NODE, h.get_node_label(*node_1), h.get_edge_label(*edge_1), h.get_node_label(node_2), *node_1, *edge_1, node_2);
			if (not is_substruct_in_g.at(node_edge_node_substruct)) {
				unmatched_substructs_.push_back(node_edge_node_substruct);
				is_deleted_node.at(*node_1) = true;
				is_deleted_edge.at(*edge_1) = true;
				is_deleted_node.at(node_2) = true;
				break;
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Partition<UserNodeLabel, UserEdgeLabel>::
contains_node_edge_node_substruct_(const GEDGraph & g, const Substruct_ & substruct) const {
	std::vector<GEDGraph::NodeID> candidate_nodes;
	for (auto node_1 = g.nodes().first; node_1 != g.nodes().second; node_1++) {
		if (g.get_node_label(*node_1) == substruct.node_label_1) {
			candidate_nodes.push_back(*node_1);
		}
	}
	if (candidate_nodes.empty()) {
		return false;
	}
	for (auto node_1 = candidate_nodes.begin(); node_1 != candidate_nodes.end(); node_1++) {
		for (auto edge_1 = g.incident_edges(*node_1).first; edge_1 != g.incident_edges(*node_1).second; edge_1++) {
			if (g.get_edge_label(*edge_1) == substruct.edge_label_1) {
				if (g.get_node_label(g.head(*edge_1)) == substruct.node_label_2) {
					return true;
				}
			}
		}
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Partition<UserNodeLabel, UserEdgeLabel>::
check_node_edge_edge_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::NodeID, bool> & is_deleted_node, std::map<GEDGraph::EdgeID, bool> & is_deleted_edge) {
	for (auto node_1 = h.nodes().first; node_1 != h.nodes().second; node_1++) {
		if (is_deleted_node.at(*node_1)) {
			continue;
		}
		for (auto edge_1 = h.incident_edges(*node_1).first; edge_1 != h.incident_edges(*node_1).second; edge_1++) {
			if (is_deleted_edge.at(*edge_1)) {
				continue;
			}
			bool found_unmatched_substruct{false};
			for (auto edge_2 = h.incident_edges(*node_1).first; edge_2 != h.incident_edges(*node_1).second; edge_2++) {
				if ((*edge_1 == *edge_2) or is_deleted_edge.at(*edge_2)) {
					continue;
				}
				Substruct_ node_edge_edge_substruct(NODE_EDGE_EDGE, h.get_node_label(*node_1), h.get_edge_label(*edge_1), h.get_edge_label(*edge_2), *node_1, *edge_1, GEDGraph::dummy_node(), *edge_2);
				if (not is_substruct_in_g.at(node_edge_edge_substruct)) {
					unmatched_substructs_.push_back(node_edge_edge_substruct);
					is_deleted_node.at(*node_1) = true;
					is_deleted_edge.at(*edge_1) = true;
					is_deleted_edge.at(*edge_2) = true;
					found_unmatched_substruct = true;
					break;
				}
			}
			if (found_unmatched_substruct) {
				break;
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Partition<UserNodeLabel, UserEdgeLabel>::
contains_node_edge_edge_substruct_(const GEDGraph & g, const Substruct_ & substruct) const {
	std::vector<GEDGraph::NodeID> candidate_nodes;
	for (auto node_1 = g.nodes().first; node_1 != g.nodes().second; node_1++) {
		if (g.get_node_label(*node_1) == substruct.node_label_1) {
			candidate_nodes.push_back(*node_1);
		}
	}
	if (candidate_nodes.empty()) {
		return false;
	}
	for (auto node_1 = candidate_nodes.begin(); node_1 != candidate_nodes.end(); node_1++) {
		for (auto edge_1 = g.incident_edges(*node_1).first; edge_1 != g.incident_edges(*node_1).second; edge_1++) {
			if (g.get_edge_label(*edge_1) == substruct.edge_label_1) {
				for (auto edge_2 = g.incident_edges(*node_1).first; edge_2 != g.incident_edges(*node_1).second; edge_2++) {
					if (*edge_1 == *edge_2) {
						continue;
					}
					if (g.get_edge_label(*edge_2) == substruct.edge_label_2) {
						return true;
					}
				}
				break;
			}
		}
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
Partition<UserNodeLabel, UserEdgeLabel>::
Substruct_ ::
Substruct_(SubstructType_ type, LabelID label_1, LabelID label_2, LabelID label_3, GEDGraph::NodeID node_1, GEDGraph::EdgeID edge_1, GEDGraph::NodeID node_2, GEDGraph::EdgeID edge_2) :
type{type},
node_label_1{dummy_label()},
node_label_2{dummy_label()},
edge_label_1{dummy_label()},
edge_label_2{dummy_label()},
node_1{node_1},
node_2{node_2},
edge_1{edge_1},
edge_2{edge_2}{
	switch (type) {
	case NODE:
		node_label_1 = label_1;
		break;
	case EDGE:
		edge_label_1 = label_1;
		break;
	case NODE_EDGE:
		node_label_1 = label_1;
		edge_label_1 = label_2;
		break;
	case NODE_EDGE_NODE:
		node_label_1 = label_1;
		edge_label_1 = label_2;
		node_label_2 = label_3;
		break;
	case NODE_EDGE_EDGE:
		node_label_1 = label_1;
		edge_label_1 = label_2;
		edge_label_2 = label_3;
		break;
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Partition<UserNodeLabel, UserEdgeLabel>::
Substruct_ ::
operator<(const Substruct_ & rhs) const {
	if (node_label_1 < rhs.node_label_1) {
		return true;
	}
	if (node_label_1 > rhs.node_label_1) {
		return false;
	}
	if (edge_label_1 < rhs.edge_label_1) {
		return true;
	}
	if (edge_label_1 > rhs.edge_label_1) {
		return false;
	}
	if (node_label_2 < rhs.node_label_2) {
		return true;
	}
	if (node_label_2 > rhs.node_label_2) {
		return false;
	}
	if (edge_label_2 < rhs.edge_label_2) {
		return true;
	}
	if (edge_label_2 > rhs.edge_label_2) {
		return false;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
const std::vector<typename Partition<UserNodeLabel, UserEdgeLabel>::Substruct_> &
Partition<UserNodeLabel, UserEdgeLabel>::
get_unmatched_substructs_() const {
	return unmatched_substructs_;
}

}

#endif /* SRC_METHODS_PARTITION_IPP_ */
