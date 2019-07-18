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
 * @file edit_costs.hpp
 * @brief ged::EditCosts class declaration.
 */

#ifndef SRC_EDIT_COSTS_EDIT_COSTS_HPP_
#define SRC_EDIT_COSTS_EDIT_COSTS_HPP_

#include "../env/common_types.hpp"

namespace ged {

/*!
 * @brief Abstract class for defining edit cost functions.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class EditCosts {

public:

	/*!
	 * @brief Virtual destructor.
	 * @note Must be implemented by derived classes of ged::EditCosts.
	 */
	virtual ~EditCosts() = 0;

	/*!
	 * @brief Constructor.
	 */
	EditCosts();

	/*!
	 * @brief Node insertions cost function.
	 * @param[in] node_label A node label.
	 * @return The cost of inserting a node with label @p node_label.
	 * @note Must be implemented by derived classes of ged::EditCosts.
	 */
	virtual double node_ins_cost_fun(const UserNodeLabel & node_label) const;

	/*!
	 * @brief Node deletion cost function.
	 * @param[in] node_label A node label.
	 * @return The cost of deleting a node with label @p node_label.
	 * @note Must be implemented by derived classes of ged::EditCosts.
	 */
	virtual double node_del_cost_fun(const UserNodeLabel & node_label) const;

	/*!
	 * @brief Node relabeling cost function.
	 * @param[in] node_label_1 A node label.
	 * @param[in] node_label_2 A node label.
	 * @return The cost of changing a node's label from @p node_label_1 to @p node_label_2.
	 * @note Must be implemented by derived classes of ged::EditCosts.
	 */
	virtual double node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const;

	/*!
	 * @brief Computes a node label's representation as a real-valued vector.
	 * @param[in] node_label A node label.
	 * @param[out] vector_representation The node label's vector representation.
	 * @note May be implemented by derived classes of ged::EditCosts.
	 */
	virtual void vectorize_node_label(const UserNodeLabel & node_label, std::vector<double> & vector_representation) const;

	/*!
	 * @brief Computes median node label.
	 * @param[in] node_labels The node labels whose median should be computed.
	 * @return Median of the given node labels.
	 * @note Must be implemented if the nodes are labeled and ged::MedianGraphEstimator should be used with the edit costs.
	 */
	virtual UserNodeLabel median_node_label(const std::vector<UserNodeLabel> & node_labels) const;

	/*!
	 * @brief Edge insertion cost function.
	 * @param[in] edge_label An edge label.
	 * @return The cost of inserting an edge with label @p edge_label.
	 * @note Must be implemented by derived classes of ged::EditCosts.
	 */
	virtual double edge_ins_cost_fun(const UserEdgeLabel & edge_label) const;

	/*!
	 * @brief Edge deletion cost function.
	 * @param[in] edge_label An edge label.
	 * @return The cost of deleting an edge with label @p edge_label.
	 * @note Must be implemented by derived classes of ged::EditCosts.
	 */
	virtual double edge_del_cost_fun(const UserEdgeLabel & edge_label) const;

	/*!
	 * @brief Edge relabeling cost function.
	 * @param[in] edge_label_1 An edge label.
	 * @param[in] edge_label_2 An edge label.
	 * @return The cost of changing an edge's label from @p edge_label_1 to @p edge_label_2.
	 * @note Must be implemented by derived classes of ged::EditCosts.
	 */
	virtual double edge_rel_cost_fun(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const;

	/*!
	 * @brief Computes an edge label's representation as a real-valued vector.
	 * @param[in] edge_label An edge label.
	 * @param[out] vector_representation The edge label's vector representation.
	 * @note May be implemented by derived classes of ged::EditCosts.
	 */
	virtual void vectorize_edge_label(const UserEdgeLabel & edge_label, std::vector<double> & vector_representation) const;

	/*!
	 * @brief Computes median edge label.
	 * @param[in] edge_labels The edge labels whose median should be computed.
	 * @return Median of the given edge labels.
	 * @note Must be implemented if the edges are labeled and ged::MedianGraphEstimator should be used with the edit costs.
	 */
	virtual UserEdgeLabel median_edge_label(const std::vector<UserEdgeLabel> & edge_labels) const;
};

}

#include "edit_costs.ipp"

#endif /* SRC_EDIT_COSTS_EDIT_COSTS_HPP_ */
