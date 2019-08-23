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
 * @file constant.hpp
 * @brief ged::Constant class declaration.
 */

#ifndef SRC_EDIT_COSTS_CONSTANT_HPP_
#define SRC_EDIT_COSTS_CONSTANT_HPP_

#include "edit_costs.hpp"

namespace ged {

/*!
 * @brief Implements constant edit cost functions.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Constant : public EditCosts<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Constant();

	/*!
	 * @brief Constructor.
	 * @param[in] node_ins_cost
	 * @param[in] node_del_cost
	 * @param[in] node_rel_cost
	 * @param[in] edge_ins_cost
	 * @param[in] edge_del_cost
	 * @param[in] edge_rel_cost
	 * @note Calling the constructor with the default arguments constructs uniform edit costs.
	 */
	Constant(double node_ins_cost = 1, double node_del_cost = 1, double node_rel_cost = 1, double edge_ins_cost = 1, double edge_del_cost = 1, double edge_rel_cost = 1);

	virtual double node_ins_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_del_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const final;

	virtual double edge_ins_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_del_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_rel_cost_fun(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const final;

private:

	double node_ins_cost_;

	double node_del_cost_;

	double node_rel_cost_;

	double edge_ins_cost_;

	double edge_del_cost_;

	double edge_rel_cost_;
};

}

#include "constant.ipp"

#endif /* SRC_EDIT_COSTS_CONSTANT_HPP_ */
