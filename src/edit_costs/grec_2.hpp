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
 * @file grec_2.hpp
 * @brief ged::GREC2 class declaration.
 */

#ifndef SRC_EDIT_COSTS_GREC_2_HPP_
#define SRC_EDIT_COSTS_GREC_2_HPP_

#include "../util/lsape_solver.hpp"
#include "edit_costs.hpp"

namespace ged {

/*!
 * @brief Edit cost functions for the dataset GREC.
 * @details The graphs contained in the GREC dataset represent line drawings of electronic or architectural symbols.
 * The contained graphs have around 20 vertices that are attributed with integer-valued Euclidian coordinates
 * (named "x" and "y"), as well as with a string type (named "type"). Edges are attributed with an integer named "frequency"
 * that takes the values 1 and 2, strings named "type0" and "angle0", and strings named "type1" and "angle1" if the frequency equals 2.
 * The GREC dataset is contained in the IAM graph database repository which can be downloaded from http://www.fki.inf.unibe.ch/databases/iam-graph-database:
 * - K. Riesen, H. Bunke:
 *   &ldquo;IAM graph database repository for graph based pattern recognition and machine learning&rdquo;,
 *   https://doi.org/10.1007/978-3-540-89689-0_33
 *
 * Implements the edit costs suggested in:
 * - Z. Abu-Aisheh, R. Raveaux, J.-Y. Ramel:
 *   &ldquo;A graph database repository and performance evaluation metrics for graph edit distance&rdquo;,
 *   https://doi.org/10.1007/978-3-319-18224-7_14
 * - K. Riesen, H. Bunke:
 *   &ldquo;Graph data&rdquo;, in: *Graph Classification and Clustering Based on Vector Space Embedding*,
 *   https://doi.org/10.1142/9789814304726_0004
 */
template<class UserNodeLabel, class UserEdgeLabel>
class GREC2 : public EditCosts<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~GREC2();

	/*!
	 * @brief Constructor.
	 * @param[in] node_ins_del_cost Cost for deleting or inserting nodes.
	 * @param[in] edge_ins_del_cost Cost for deleting or inserting edges.
	 * @param[in] alpha Importance of node edit operations vs. importance of edge edit operations.
	 * @note Calling the constructor with the default arguments constructs the edit costs suggested in https://doi.org/10.1016/j.patrec.2017.10.007.
	 */
	GREC2(double node_ins_del_cost = 90, double edge_ins_del_cost = 15, double alpha = 0.5);

	virtual double node_ins_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_del_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const final;

	virtual double edge_ins_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_del_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_rel_cost_fun(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const final;

private:

	double node_ins_del_cost_;

	double edge_ins_del_cost_;

	double alpha_;
};

}

#include "grec_2.ipp"

#endif /* SRC_EDIT_COSTS_GREC_2_HPP_ */
