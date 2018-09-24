/*!
 * @file node_map.ipp
 * @brief ged::NodeMap class definition.
 */

#ifndef SRC_ENV_NODE_MAP_IPP_
#define SRC_ENV_NODE_MAP_IPP_

namespace ged {

NodeMap::
NodeMap() :
forward_map_(),
backward_map_(),
induced_cost_{std::numeric_limits<double>::infinity()} {}

NodeMap::
NodeMap(const NodeMap & matching) :
forward_map_(matching.forward_map_),
backward_map_(matching.backward_map_),
induced_cost_(matching.induced_cost_) {}

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
	forward_map_.clear();
	backward_map_.clear();
}

bool
NodeMap::
empty() const {
	return forward_map_.empty() and backward_map_.empty();
}

bool
NodeMap::
complete(const GEDGraph & g, const GEDGraph & h) const {
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		if (forward_map_.find(*node) == forward_map_.end()) {
			return false;
		}
	}
	for (auto node = h.nodes().first; node != h.nodes().second; node++) {
		if (backward_map_.find(*node) == backward_map_.end()) {
			return false;
		}
	}
	return true;
}

GEDGraph::NodeID
NodeMap::
image(GEDGraph::NodeID node) const {
	try {
		return forward_map_.at(node);
	}
	catch (const std::out_of_range& oor) {
		throw Error("The node map does not contain a substitution or a deletion for node " + std::to_string(node) + ".");
	}
	return GEDGraph::dummy_node();
}

void
NodeMap::
erase_image(GEDGraph::NodeID node) {
	forward_map_.erase(node);
}

GEDGraph::NodeID
NodeMap::
pre_image(GEDGraph::NodeID node) const {
	try {
		return backward_map_.at(node);
	}
	catch (const std::out_of_range& oor) {
		throw Error("The node map does not contain a substitution or an insertion for node " + std::to_string(node) + ".");
	}
	return GEDGraph::dummy_node();
}

void
NodeMap::
erase_pre_image(GEDGraph::NodeID node) {
	backward_map_.erase(node);
}

void
NodeMap::
as_relation(std::vector<Assignment> & relation) const {
	relation.clear();
	for (auto assignment = forward_map_.begin(); assignment != forward_map_.end(); assignment++) {
		relation.emplace_back(assignment->first, assignment->second);
	}
	for (auto assignment = backward_map_.begin(); assignment != backward_map_.end(); assignment++) {
		if (assignment->second == GEDGraph::dummy_node()) {
			relation.emplace_back(assignment->second, assignment->first);
		}
	}
}

void
NodeMap::
add_assignment(GEDGraph::NodeID i, GEDGraph::NodeID k) {
	if (i != GEDGraph::dummy_node()) {
		forward_map_[i] = k;
	}
	if (k != GEDGraph::dummy_node()) {
		backward_map_[k] = i;
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
