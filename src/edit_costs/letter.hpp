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
 * @file letter.hpp
 * @brief ged::Letter class declaration.
 */

#ifndef SRC_EDIT_COSTS_LETTER_HPP_
#define SRC_EDIT_COSTS_LETTER_HPP_

#include "edit_costs.hpp"

namespace ged {

/*!
 * @brief Edit costs for graphs contained in Letter datasets.
 * @details The graphs contained in the Letter datasets represent the capital letters A, E, F, H, I, K, L, M, N, T, V, W, X, Y, and Z which have been distorted by three different degrees (low, medium, high).
 * Nodes are attributed with Euclidean coordinates (named "x" and "y").
 * Edges are have no attributes.
 * The Letter datasets are contained in the IAM graph database repository which can be downloaded from http://www.fki.inf.unibe.ch/databases/iam-graph-database:
 * - K. Riesen, H. Bunke:
 *   &ldquo;IAM graph database repository for graph based pattern recognition and machine learning&rdquo;,
 *   https://doi.org/10.1007/978-3-540-89689-0_33
 *
 * Implements the edit costs suggested in:
 * - K. Riesen, H. Bunke:
 *   &ldquo;Graph data&rdquo;, in: *Graph Classification and Clustering Based on Vector Space Embedding*,
 *   https://doi.org/10.1142/9789814304726_0004
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Letter : public EditCosts<UserNodeLabel, UserEdgeLabel> {
public:

	virtual ~Letter();

	/*!
	 * @brief Constructor.
	 * @param[in] euclidean_exponent Raise Euclidean node substitution cost to the power of @p exponent.
	 * @param[in] node_ins_del_cost Cost for deleting or inserting nodes.
	 * @param[in] edge_ins_del_cost Cost for deleting or inserting edges.
	 * @param[in] alpha Importance of node edit operations vs. importance of edge edit operations.
	 * @note Calling the constructor with the default arguments constructs the edit costs for Letter high suggested in https://doi.org/10.1142/9789814304726_0004.
	 * For Letter medium, the suggested arguments are <tt>node_ins_del_cost = 0.7</tt>, <tt>edge_ins_del_cost = 1.9</tt>, and <tt>alpha = 0.75</tt>.
	 * For Letter low, the suggested arguments are <tt>node_ins_del_cost = 0.3</tt>, <tt>edge_ins_del_cost = 0.1</tt>, and <tt>alpha = 0.25</tt>.
	 */
	Letter(double euclidean_exponent = 1, double node_ins_del_cost = 0.9, double edge_ins_del_cost = 1.7, double alpha = 0.75);

	virtual double node_ins_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_del_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const final;

	virtual UserNodeLabel median_node_label(const std::vector<UserNodeLabel> & node_labels) const final;

	virtual double edge_ins_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_del_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_rel_cost_fun(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const final;

private:

	double euclidean_exponent_;

	double node_ins_del_cost_;

	double edge_ins_del_cost_;

	double alpha_;
};

}

#include "letter.ipp"

#endif /* SRC_EDIT_COSTS_LETTER_HPP_ */
