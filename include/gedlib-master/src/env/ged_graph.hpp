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
 * @file  ged_graph.hpp
 * @brief ged::GEDGraph class declaration.
 */

#ifndef SRC_ENV_GED_GRAPH_HPP_
#define SRC_ENV_GED_GRAPH_HPP_

#include "common_types.hpp"
#include "matrix.hpp"

//#include "ged_graph.fwd.hpp"

namespace ged {

namespace detail {

namespace al_graph_prop {

struct node;

struct edge;

struct graph;

}

namespace am_graph_prop {

struct node;

struct edge;

struct graph;

}

typedef boost::adjacency_list<
		boost::vecS,
		boost::vecS,
		boost::undirectedS,
		al_graph_prop::node,
		al_graph_prop::edge
		> GedGraphAL;

typedef boost::adjacency_matrix<
		boost::undirectedS,
		am_graph_prop::node,
		am_graph_prop::edge
		> GedGraphAM;

namespace al_graph_prop {

struct node {
	LabelID label{invalid_label()};
};

struct edge {
	LabelID label{invalid_label()};
};

struct graph {};

}

namespace am_graph_prop {

struct node {};

struct edge {
	GedGraphAL::edge_descriptor cref;
};

struct graph {};

}

}

/*!
 * @brief The normalized input graphs used by GEDLIB. All labels are integers.
 */
class GEDGraph {

public:

	typedef std::size_t NodeID; //!< Internally used vertex ID type.

	typedef detail::GedGraphAL::edge_descriptor EdgeID; //!< Internally used edge ID type.

	typedef std::vector<GEDGraph>::size_type GraphID; //!< Type of internally used graph IDs.

	typedef std::map<NodeID, NodeID> NodeNodeMap; //!< Map that assigns node IDs to node IDs.

	typedef std::map<NodeID, std::size_t> NodeSizeTMap; //!< Map that assigns node IDs to integers.

	typedef std::map<std::size_t, NodeID> SizeTNodeMap; //!< Map that assigns node IDs to integers.

	typedef detail::GedGraphAL::out_edge_iterator incident_edge_iterator; //!< Iterator type for iterating over all the incident edges of a node in the graph.

	typedef detail::GedGraphAL::edge_iterator edge_iterator; //!< Iterator type for iterating over all edges of the graph.

	typedef detail::GedGraphAL::vertex_iterator node_iterator; //!< Iterator type for iterating over all nodes in the graph.

	/*!
	 * @brief Constructs an empty graph without meaningful ID.
	 */
	GEDGraph();

	/*!
	 * @brief Constructs an empty graph.
	 * @param[in] id The ID of the new graph.
	 */
	GEDGraph(GraphID id);

	/*!
	 * @brief Copy constructor.
	 * @param[in] ged_graph Graph that is to be copied.
	 */
	GEDGraph(const GEDGraph & ged_graph);

	/*!
	 * @brief Returns the ID of the graph.
	 * @return ID of the graph.
	 */
	GraphID id() const;

	/*!
	 * @brief Clears the graph.
	 */
	void clear();

	/*!
	 * @brief Returns a dummy node.
	 * @return ID of dummy node.
	 */
	static NodeID dummy_node();

	/*!
	 * @brief Returns an undefined node.
	 * @return ID of undefined node.
	 */
	static NodeID undefined_node();

	/*!
	 * @brief Returns a dummy edge.
	 * @return ID of dummy edge.
	 */
	static EdgeID dummy_edge();


	/*!
	 * @brief Adds a new vertex to the graph.
	 * @return The ID of the newly added vertex.
	 */
	NodeID add_node();

	/*!
	 * @brief Adds a new edge to the graph.
	 * @param[in] tail The tail of the newly added edge.
	 * @param[in] head The head of the newly added edge
	 * @return The ID of the newly added edge.
	 */
	EdgeID add_edge(NodeID tail, NodeID head);

	/*!
	 * @brief Sets a node's label.
	 * @param[in] node The ID of the node whose label should be set.
	 * @param[in] label The ID of the label that should be assigned to the node @p node.
	 */
	void set_label(NodeID node, LabelID label);

	/*!
	 * @brief Sets a edges's label.
	 * @param[in] edge The ID of the edge whose label should be set.
	 * @param[in] label The ID of the label that should be assigned to the edge @p edge.
	 */
	void set_label(EdgeID edge, LabelID label);

	/*!
	 * @brief Returns the label of a given node.
	 * @param[in] node The ID of the node whose label is to be returned.
	 * @return The ID of the label of node @p node.
	 */
	LabelID get_node_label(NodeID node) const;

	/*!
	 * @brief Returns the label of a given edge.
	 * @param[in] edge The ID of the edge whose label is to be returned.
	 * @return The ID of the label of edge @p edge.
	 */
	LabelID get_edge_label(EdgeID edge) const;

	/*!
	 * @brief Retrieves the label of an edge between two nodes or the dummy label if no edge exists.
	 * @param[in] node_1 ID of first node.
	 * @param[in] node_2 ID of second node.
	 * @return The label of the edge that connects @p node_1 with @p node_2 or the dummy label if no edge exists.
	 */
	LabelID get_edge_label(NodeID node_1, NodeID node_2) const;

	/*!
	 * @brief Prepares the adjacency matrix to allow quick retrieval of edges.
	 * @warning Calls to the member add_edge() invalidate the matrix.
	 */
	void setup_adjacency_matrix();

	/*!
	 * @brief Returns the tail of an edge.
	 * @param[in] edge The ID of the edge whose tail is to be returned.
	 * @return The ID of the tail of edge's internal representation as directed edge.
	 * @note Since GEDGraphs are undirected, the tail of an edge is not well defined.
	 * However, given the same edge ID, tail() always returns the same node.
	 */
	NodeID tail(EdgeID edge) const;

	/*!
	 * @brief Returns the head of an edge.
	 * @param[in] edge The ID of the edge whose head is to be returned.
	 * @return The ID of the head of edge's internal representation as directed edge.
	 * @note Since GEDGraphs are undirected, the head of an edge is not well defined.
	 * However, given the same edge ID, head() always returns the same node.
	 */
	NodeID head(EdgeID edge) const;


	/*!
	 * @brief Retrieves an edge from its incident nodes.
	 * @param[in] tail First incident node of the edge.
	 * @param[in] head Second incident node of the edge.
	 * @return The ID of the edge that connects @p tail with @p head, if the edge exists. Otherwise, dummy_edge() is returned.
	 * @warning This function uses the adjacency matrix. Ensure you called setup_adjacency_matrix() before calling.
	 */
	EdgeID get_edge(NodeID tail, NodeID head) const;


	/*!
	 * @brief Retrieves an edge from its incident nodes.
	 * @param[in] tail First incident node of the edge.
	 * @param[in] head Second incident node of the edge.
	 * @return The ID of the edge that connects @p tail with @p head, if the edge exists. Otherwise, dummy_edge() is returned.
	 * @note This function uses the adjacency list. It is hence slower than get_edge() but does not expect that setup_adjacency_matrix() has been called.
	 */
	EdgeID safe_get_edge(NodeID tail, NodeID head) const;


	/*!
	 * @brief Checks if an edge exists
	 * @param[in] tail First node.
	 * @param[in] head Second node.
	 * @return Boolean @p true if there is an edge that connects @p tail with @p head, and @p false otherwise.
	 * @warning This function uses the adjacency matrix. Ensure you called setup_adjacency_matrix() before calling.
	 */
	bool is_edge(NodeID tail, NodeID head) const;

	/*!
	 * @brief Checks if an edge exists
	 * @param[in] tail First node.
	 * @param[in] head Second node.
	 * @return Boolean @p true if there is an edge that connects @p tail with @p head, and @p false otherwise.
	 * @note This function uses the adjacency list. It is hence slower than is_edge() but does not expect that setup_adjacency_matrix() has been called.
	 */
	bool safe_is_edge(NodeID tail, NodeID head) const;

	/*!
	 * @brief Returns node degree.
	 * @param[in] node Input node.
	 * @return Degree of node @p node.
	 */
	std::size_t degree(NodeID node) const;

	/*!
	 * @brief Returns the number of nodes.
	 * @return Number of nodes contained in the graph.
	 */
	std::size_t num_nodes() const;

	/*!
	 * @brief Returns the number of edges.
	 * @return Number of edges contained in the graph.
	 */
	std::size_t num_edges() const;

	/*!
	 * @brief Checks if graph has already been initialized.
	 * @return Boolean @p true if the adjacency matrix is up to date, @p false otherwise.
	 */
	bool initialized() const;

	/*!
	 * @brief Declare graph as un-initialized.
	 */
	void un_init();

	/*!
	 * \brief Returns the maximum degree of the graph.
	 */
	std::size_t maxdeg() const;

	/*!
	 * @brief Provides access to all nodes in the graph.
	 * @return A pair of iterators whose first member points to the first node in the graph and whose second member points past the last node in the graph.
	 */
	std::pair<node_iterator, node_iterator> nodes() const;

	/*!
	 * @brief Provides access to all edge in the graph.
	 * @return A pair of iterators whose first member points to the first edge in the graph and whose second member points past the last edge in the graph.
	 */
	std::pair<edge_iterator, edge_iterator> edges() const;

	/*!
	 * @brief Provides access to all incident edges of a node.
	 * @param[in] node ID of the node whose incident edges should be returned.
	 * @return A pair of iterators whose first member points to the first incident edge of the node @p node and whose second member points past the last incident edge of the node @p node.
	 */
	std::pair<incident_edge_iterator, incident_edge_iterator> incident_edges(NodeID node) const;

	/*!
	 * @brief Returns the internal adjacency list representation.
	 * @return Internal adjacency list representation.
	 */
	const detail::GedGraphAL & get_adjacency_list() const;

private:

	detail::GedGraphAL alist_;

	detail::GedGraphAM amatrix_;

	GraphID id_;

	bool initialized_;
};

/*!
 * @brief Streams a graph.
 * @param[in,out] os Output stream.
 * @param[in] graph The graph that should be streamed.
 * @return The output stream @p os.
 * @relates ged::GEDGraph
 */
std::ostream & operator<<(std::ostream & os, const GEDGraph & graph);

}

#include "ged_graph.ipp"

#endif /* SRC_ENV_GED_GRAPH_HPP_ */
