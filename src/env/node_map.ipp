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
 * @file node_map.ipp
 * @brief ged::NodeMap class definition.
 */

#ifndef SRC_ENV_NODE_MAP_IPP_
#define SRC_ENV_NODE_MAP_IPP_

namespace ged {

NodeMap::
NodeMap(std::size_t num_nodes_g, std::size_t num_nodes_h) :
forward_map_(num_nodes_g, GEDGraph::undefined_node()),
backward_map_(num_nodes_h, GEDGraph::undefined_node()),
induced_cost_{std::numeric_limits<double>::infinity()} {}

NodeMap::
NodeMap(const NodeMap & node_map) :
forward_map_(node_map.forward_map_),
backward_map_(node_map.backward_map_),
induced_cost_(node_map.induced_cost_) {}

void
NodeMap::
operator=(const NodeMap & node_map) {
	forward_map_ = node_map.forward_map_;
	backward_map_ = node_map.backward_map_;
	induced_cost_ = node_map.induced_cost_;
}

void
NodeMap::
clear() {
	for (auto & node_id : forward_map_) {
		node_id = GEDGraph::undefined_node();
	}
	for (auto & node_id : backward_map_) {
		node_id = GEDGraph::undefined_node();
	}
}

bool
NodeMap::
empty() const {
	for (const auto & node_id : forward_map_) {
		if (node_id != GEDGraph::undefined_node()) {
			return false;
		}
	}
	for (const auto & node_id : backward_map_) {
		if (node_id != GEDGraph::undefined_node()) {
			return false;
		}
	}
	return true;
}

std::size_t
NodeMap::
num_source_nodes() const {
	return forward_map_.size();
}

std::size_t
NodeMap::
num_target_nodes() const {
	return backward_map_.size();
}

bool
NodeMap::
complete(const GEDGraph & g, const GEDGraph & h) const {
	for (const auto & node_id : forward_map_) {
		if (node_id == GEDGraph::undefined_node()) {
			return false;
		}
	}
	for (const auto & node_id : backward_map_) {
		if (node_id == GEDGraph::undefined_node()) {
			return false;
		}
	}
	return true;
}


GEDGraph::NodeID
NodeMap::
image(GEDGraph::NodeID node) const {
	if (node < forward_map_.size()) {
		return forward_map_.at(node);
	}
	else {
		throw Error("The node with ID " + std::to_string(node) + " is not contained in the source nodes of the node map.");
	}
	return GEDGraph::undefined_node();
}

void
NodeMap::
erase_image(GEDGraph::NodeID node) {
	if (node < forward_map_.size()) {
		GEDGraph::NodeID image{forward_map_.at(node)};
		forward_map_[node] = GEDGraph::undefined_node();
		if (image < backward_map_.size()) {
			backward_map_[image] = GEDGraph::undefined_node();
		}
	}
	else {
		throw Error("The node with ID " + std::to_string(node) + " is not contained in the source nodes of the node map.");
	}
}

GEDGraph::NodeID
NodeMap::
pre_image(GEDGraph::NodeID node) const {
	if (node < backward_map_.size()) {
		return backward_map_.at(node);
	}
	else {
		throw Error("The node with ID " + std::to_string(node) + " is not contained in the target nodes of the node map.");
	}
	return GEDGraph::undefined_node();
}

std::vector<GEDGraph::NodeID>
NodeMap::
get_forward_map() const {
	return forward_map_;
}

std::vector<GEDGraph::NodeID>
NodeMap::
get_backward_map() const {
	return backward_map_;
}

void
NodeMap::
erase_pre_image(GEDGraph::NodeID node) {
	if (node < backward_map_.size()) {
		GEDGraph::NodeID pre_image{backward_map_.at(node)};
		backward_map_[node] = GEDGraph::undefined_node();
		if (pre_image < forward_map_.size()) {
			forward_map_[pre_image] = GEDGraph::undefined_node();
		}
	}
	else {
		throw Error("The node with ID " + std::to_string(node) + " is not contained in the target nodes of the node map.");
	}
}

void
NodeMap::
as_relation(std::vector<Assignment> & relation) const {
	relation.clear();
	GEDGraph::NodeID i, k;
	for (i = 0; i < forward_map_.size(); i++) {
		k = forward_map_.at(i);
		if (k != GEDGraph::undefined_node()) {
			relation.emplace_back(i, k);
		}
	}
	for (k = 0; k < backward_map_.size(); k++) {
		i = backward_map_.at(k);
		if (i == GEDGraph::dummy_node()) {
			relation.emplace_back(i, k);
		}
	}
}

void
NodeMap::
add_assignment(GEDGraph::NodeID i, GEDGraph::NodeID k) {
	if (i != GEDGraph::dummy_node()) {
		if (i < forward_map_.size()) {
			forward_map_[i] = k;
		}
		else {
			throw Error("The node with ID " + std::to_string(i) + " is not contained in the source nodes of the node map.");
		}
	}
	if (k != GEDGraph::dummy_node()) {
		if (k < backward_map_.size()) {
			backward_map_[k] = i;
		}
		else {
			throw Error("The node with ID " + std::to_string(k) + " is not contained in the target nodes of the node map.");
		}
	}
}

void
NodeMap::
set_induced_cost(double induced_cost) {
	induced_cost_= induced_cost;
}

double
NodeMap::
induced_cost() const {
	return induced_cost_;
}

bool
NodeMap::
operator<(const NodeMap & rhs) const {
	return (induced_cost_ < rhs.induced_cost_);
}

bool
NodeMap::
operator==(const NodeMap & node_map) const {
	return ((forward_map_ == node_map.forward_map_) and (backward_map_ == node_map.backward_map_));
}

std::ostream & operator<<(std::ostream & os, const NodeMap & node_map) {
	std::vector<NodeMap::Assignment> relation;
	node_map.as_relation(relation);
	os << "{ ";
	for (const auto & assignment : relation) {
		os << assignment << " ";
	}
	os << "}";
	return os;
}

std::ostream & operator<<(std::ostream & os, const NodeMap::Assignment & assignment) {
	os << "(";
	if (assignment.first == GEDGraph::dummy_node()) {
		os << "DUMMY,";
	}
	else {
		os << assignment.first << ",";
	}
	if (assignment.second == GEDGraph::dummy_node()) {
		os << "DUMMY";
	}
	else {
		os << assignment.second << "";
	}
	os << ")";
	return os;
}

}

#endif /* SRC_ENV_NODE_MAP_IPP_ */
