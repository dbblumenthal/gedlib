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
 * @file node_map.hpp
 * @brief ged::NodeMap class declaration.
 */

#ifndef SRC_ENV_NODE_MAP_HPP_
#define SRC_ENV_NODE_MAP_HPP_

#include "ged_graph.hpp"
#include "common_types.hpp"

namespace ged {

/*!
 * @brief A class for node maps.
 * @details Node maps are relations between the node sets of two graphs augmented with one dummy node each such that each real node is assigned at most once.
 * The exact definition can be found in:
 * D. B. Blumenthal, J. Gamper.
 * On the exact computation of the graph edit distance.
 * Pattern Recogn. Lett. 2018, in press.
 */
class NodeMap {

public:

	typedef std::pair<GEDGraph::NodeID, GEDGraph::NodeID> Assignment;

	/*!
	 * @brief Default constructor.
	 * @param[in] num_nodes_g Number of nodes in first input graph.
	 * @param[in] num_nodes_h Number of nodes in second input graph.
	 */
	NodeMap(std::size_t num_nodes_g, std::size_t num_nodes_h);

	/*!
	 * @brief Copy constructor.
	 * @param[in] node_map Node map that should be copied.
	 */
	NodeMap(const NodeMap & node_map);

	/*!
	 * @brief Assignment operator.
	 * @param[in] node_map Node map that should be assigned.
	 */
	void operator=(const NodeMap & node_map);

	/*!
	 * @brief Clears the node map.
	 */
	void clear();

	/*!
	 * @brief Checks if the node map is empty.
	 * @return Boolean @p true if the map is empty and @p false otherwise.
	 */
	bool empty() const;

	/*!
	 * @brief Returns number of source nodes contained in the node map.
	 * @return Number of source nodes.
	 */
	std::size_t num_source_nodes() const;

	/*!
	 * @brief Returns number of target nodes contained in the node map.
	 * @return Number of target nodes.
	 */
	std::size_t num_target_nodes() const;

	/*!
	 * @brief Checks if the node map is complete.
	 * @param[in] g Source graph.
	 * @param[in] h Target graph.
	 * @return Boolean @p true if and only if, for all nodes in the graphs @p g and @p h, the node map contains a substitution, an insertion, or a deletion.
	 */
	bool complete(const GEDGraph & g, const GEDGraph & h) const;

	/*!
	 * @brief Returns image of a node.
	 * @param[in] node Node whose image is to be returned.
	 * @return Node to which node @p node is assigned.
	 */
	GEDGraph::NodeID image(GEDGraph::NodeID node) const;

	/*!
	 * @brief Returns pre-image of a node.
	 * @param[in] node Node whose pre_image is to be returned.
	 * @return Node which is assigned to node @p node.
	 */
	GEDGraph::NodeID pre_image(GEDGraph::NodeID node) const;

	/*!
	 * @brief Returns the forward map.
	 * @return Forward map.
	 */
	std::vector<GEDGraph::NodeID> get_forward_map() const;

	/*!
	 * @brief Returns the backward map.
	 * @return Backward map.
	 */
	std::vector<GEDGraph::NodeID> get_backward_map() const;

	/*!
	 * @brief Constructs the representation as relation.
	 * @param[out] relation Contains the pairs of nodes in the relation.
	 * The first member of a pair contained in the relation is a node in the first graph or the dummy node,
	 * the second member is a node in the second graph or the dummy node.
	 */
	void as_relation(std::vector<Assignment> & relation) const;

	/*!
	 * @brief Returns the induced cost.
	 * @return The cost induced by the node map. Equals std::numeric_limits<double>::infinity() if the induced has not been set by call to set_induced_cost().
	 */
	double induced_cost() const;

	/*!
	 * @brief Add node substitution, insertion, or deletion to the node map.
	 * @param[in] i Source node. If ged::GEDGraph::dummy_node(), the method adds a node insertion to the node map.
	 * @param[in] k Target node. If ged::GEDGraph::dummy_node(), the method adds a node deletion to the node map.
	 * @note If both @p i and @p k equal ged::GEDGraph::dummy_node(), calling this method has no effect.
	 */
	void add_assignment(GEDGraph::NodeID i, GEDGraph::NodeID k);

	/*!
	 * @brief Erases image of a node.
	 * @param[in] node Node whose image is to be erased.
	 */
	void erase_image(GEDGraph::NodeID node);

	/*!
	 * @brief Erases pre-image of a node.
	 * @param[in] node Node whose pre-image is to be erased.
	 */
	void erase_pre_image(GEDGraph::NodeID node);

	/*!
	 * @brief Sets the induced cost of the node map.
	 * @param[in] induced_cost The induced cost of the node map.
	 */
	void set_induced_cost(double induced_cost);

	/*!
	 * @brief Compares this node map to another node map.
	 * @param[in] node_map The node to which this node map should be compared.
	 * @return Boolean @p true if the induced cost is smaller than the induced cost of @p node_map and @p false otherwise.
	 */
	bool operator<(const NodeMap & node_map) const;

	/*!
	 * @brief Checks if two node maps are the same.
	 * @param[in] node_map The node map to which this node map should be compared.
	 * @return Boolean @p true if both node maps are the same and @p false otherwise.
	 */
	bool operator==(const NodeMap & node_map) const;

private:

	std::vector<GEDGraph::NodeID> forward_map_;

	std::vector<GEDGraph::NodeID> backward_map_;

	double induced_cost_;
};

/*!
 * @brief Streams a node map.
 * @param[in,out] os Output stream.
 * @param[in] node_map The node map that should be streamed.
 * @return The output stream @p os.
 * @relates ged::NodeMap
 */
std::ostream & operator<<(std::ostream & os, const NodeMap & node_map);

/*!
 * @brief Streams an assignment.
 * @param[in,out] os Output stream.
 * @param[in] assignment The assignment that should be streamed.
 * @return The output stream @p os.
 * @relates ged::NodeMap::Assignment
 */
std::ostream & operator<<(std::ostream & os, const NodeMap::Assignment & assignment);

}

#include "node_map.ipp"

#endif /* SRC_ENV_NODE_MAP_HPP_ */
