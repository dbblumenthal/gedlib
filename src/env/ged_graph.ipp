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
 * @file  ged_graph.ipp
 * @brief GEDGraph class definition.
 */

#ifndef SRC_ENV_GED_GRAPH_IPP_
#define SRC_ENV_GED_GRAPH_IPP_


namespace ged {

GEDGraph ::
GEDGraph():
alist_{},
amatrix_(1),
id_(undefined()),
initialized_{false} {}

GEDGraph ::
GEDGraph(GEDGraph::GraphID id):
alist_{},
amatrix_(1),
id_{id},
initialized_{false} {}

GEDGraph ::
GEDGraph(const GEDGraph & ged_graph):
alist_(ged_graph.alist_),
amatrix_(ged_graph.amatrix_),
id_{ged_graph.id_},
initialized_{ged_graph.initialized_} {}

GEDGraph :: NodeID
GEDGraph ::
dummy_node() {
	return std::numeric_limits<NodeID>::max() - 1;
}

GEDGraph :: NodeID
GEDGraph ::
undefined_node() {
	return std::numeric_limits<NodeID>::max();
}

GEDGraph :: EdgeID
GEDGraph ::
dummy_edge() {
	return std::numeric_limits<EdgeID>::max();
}

GEDGraph :: NodeID
GEDGraph ::
add_node() {
	initialized_ = false;
	NodeID new_node{boost::add_vertex(alist_)};
	if (new_node >= dummy_node()) {
		throw Error("Cannot add node. Maximal number of nodes reached.");
	}
	alist_[new_node].label = invalid_label();
	return new_node;
}

GEDGraph :: EdgeID
GEDGraph ::
add_edge(NodeID tail, NodeID head) {
	initialized_ = false;
	auto alist_ed = boost::add_edge(tail, head, alist_);
	alist_[alist_ed.first].label = invalid_label();
	return alist_ed.first;
}

void
GEDGraph ::
setup_adjacency_matrix() {
	detail::GedGraphAM tmp(num_nodes());
	for (auto eitr = boost::edges(alist_); eitr.first != eitr.second; ++eitr.first) {
		EdgeID e(*eitr.first);
		NodeID source {boost::source(e, alist_)};
		NodeID target {boost::target(e, alist_)};

		auto am_edge = boost::add_edge(source, target, tmp);
		if (am_edge.second == false) {
			throw Error("Parallel edge detected!");
		}
		tmp[am_edge.first].cref = e;
	}
	std::swap(tmp, amatrix_);
	initialized_ = true;
}

void
GEDGraph ::
set_label(NodeID v, LabelID l_id) {
	alist_[v].label = l_id;
}

void
GEDGraph ::
set_label(EdgeID e, LabelID l_id) {
	alist_[e].label = l_id;
}

LabelID
GEDGraph ::
get_node_label(NodeID node) const {
	if (node == dummy_node()) {
		return dummy_label();
	}
	return alist_[node].label;
}

LabelID
GEDGraph ::
get_edge_label(EdgeID e) const {
	return alist_[e].label;
}

LabelID
GEDGraph ::
get_edge_label(NodeID node_1, NodeID node_2) const {
	if (is_edge(node_1, node_2)) {
		return get_edge_label(get_edge(node_1, node_2));
	}
	return dummy_label();
}

std::pair<GEDGraph::incident_edge_iterator, GEDGraph::incident_edge_iterator>
GEDGraph ::
incident_edges(NodeID v) const {
	return boost::out_edges(v, alist_);
}

const detail::GedGraphAL &
GEDGraph ::
get_adjacency_list() const {
	return alist_;
}

std::size_t
GEDGraph ::
degree(NodeID v) const {
	return static_cast<std::size_t>(boost::out_degree(v, alist_));
}

std::size_t
GEDGraph ::
maxdeg() const {
	std::size_t maxdeg{0};
	for (auto node = this->nodes().first; node != this->nodes().second; node++) {
		maxdeg = std::max(this->degree(*node), maxdeg);
	}
	return static_cast<std::size_t>(maxdeg);
}

GEDGraph :: NodeID
GEDGraph ::
tail(EdgeID e) const {
	return boost::source(e, alist_);
}

GEDGraph::GraphID
GEDGraph::
id() const {
	return id_;
}

void
GEDGraph::
clear() {
	initialized_ = false;
	alist_.clear();
	detail::GedGraphAM tmp(1);
	std::swap(tmp, amatrix_);
}

GEDGraph :: NodeID
GEDGraph ::
head(EdgeID e) const {
	return boost::target(e, alist_);
}

std::pair<GEDGraph::node_iterator, GEDGraph::node_iterator>
GEDGraph ::
nodes() const {
	return boost::vertices(alist_);
}

std::size_t
GEDGraph ::
num_nodes() const {
	return static_cast<std::size_t>(boost::num_vertices(alist_));
}

std::pair<GEDGraph::edge_iterator, GEDGraph::edge_iterator>
GEDGraph ::
edges() const {
	return boost::edges(alist_);
}

std::size_t
GEDGraph ::
num_edges() const {
	return static_cast<std::size_t>(boost::num_edges(alist_));
}

bool
GEDGraph::
initialized() const {
	return initialized_;
}

void
GEDGraph::
un_init() {
	initialized_ = false;
}

GEDGraph::EdgeID
GEDGraph::
get_edge(NodeID source, NodeID target) const {
	//if (source >= num_nodes() or target >= num_nodes()) {
	//	return dummy_edge();
	//}
	auto am_edge = boost::edge(source, target, amatrix_);
	if (am_edge.second == false) {
		return dummy_edge();
	}
	return amatrix_[am_edge.first].cref;
}

GEDGraph::EdgeID
GEDGraph::
safe_get_edge(NodeID source, NodeID target) const {
	for (auto iedges = incident_edges(source); iedges.first != iedges.second; ++iedges.first) {
		EdgeID e(*iedges.first);
		NodeID n{this->head(e)};
		if (n == target) return e;
	}
	return dummy_edge();
}

bool
GEDGraph ::
is_edge(NodeID source, NodeID target) const {
	if (source >= this->num_nodes() or target >= this->num_nodes()) {
		return false;
	}
	auto am_edge = boost::edge(source, target, amatrix_);
	return am_edge.second;
}

bool
GEDGraph ::
safe_is_edge(NodeID source, NodeID target) const {
	for (auto iedges = incident_edges(source); iedges.first != iedges.second; ++iedges.first) {
		EdgeID e(*iedges.first);
		NodeID n{this->head(e)};
		if (n == target) return true;
	}
	return false;
}

std::ostream & operator<<(std::ostream & os, const GEDGraph & graph) {
	os << "ID = " << graph.id() <<  ", number of nodes = " << graph.num_nodes() << ", number of edges = " << graph.num_edges();
	os << "\nadjacency matrix:\n    ";
	for (auto head = graph.nodes().first; head != graph.nodes().second; head++) {
		os << "\033[1m" << std::setw(3) << std::right << graph.get_node_label(*head) << "\033[0m" << " ";
	}
	os << "\n";
	for (auto tail = graph.nodes().first; tail != graph.nodes().second; tail++) {
		os << "\033[1m" << std::setw(3) << std::right << graph.get_node_label(*tail) << "\033[0m" << " ";
		for (auto head = graph.nodes().first; head != graph.nodes().second; head++) {
			if (graph.initialized()) {
				os << "\033[1m" << std::setw(3) << std::right << graph.get_edge_label(*tail, *head) << "\033[0m" << " ";
			}
			else {
				if (graph.safe_is_edge(*tail, *head)) {
					os << "\033[1m" << std::setw(3) << std::right << graph.get_edge_label(graph.safe_get_edge(*tail, *head)) << "\033[0m" << " ";
				}
				else {
					os << std::setw(3) << std::right << 0 << " ";
				}
			}
		}
		os << "\n";
	}
	return os;
}

}

#endif /* SRC_ENV_GED_GRAPH_IPP_ */
