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
 * @file constant.ipp
 * @brief ged::Constant class definition.
 */

#ifndef SRC_EDIT_COSTS_CONSTANT_IPP_
#define SRC_EDIT_COSTS_CONSTANT_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
Constant<UserNodeLabel, UserEdgeLabel>::
~Constant() {}

template<class UserNodeLabel, class UserEdgeLabel>
Constant<UserNodeLabel, UserEdgeLabel>::
Constant(double node_ins_cost, double node_del_cost, double node_rel_cost, double edge_ins_cost, double edge_del_cost, double edge_rel_cost) :
node_ins_cost_{node_ins_cost},
node_del_cost_{node_del_cost},
node_rel_cost_{node_rel_cost},
edge_ins_cost_{edge_ins_cost},
edge_del_cost_{edge_del_cost},
edge_rel_cost_{edge_rel_cost} {}

template<class UserNodeLabel, class UserEdgeLabel>
double
Constant<UserNodeLabel, UserEdgeLabel>::
node_ins_cost_fun(const UserNodeLabel & node_label) const {
	return node_ins_cost_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Constant<UserNodeLabel, UserEdgeLabel>::
node_del_cost_fun(const UserNodeLabel & node_label) const {
	return node_del_cost_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Constant<UserNodeLabel, UserEdgeLabel>::
node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const {
	if (node_label_1 != node_label_2) {
		return node_rel_cost_;
	}
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Constant<UserNodeLabel, UserEdgeLabel>::
edge_ins_cost_fun(const UserEdgeLabel & edge_label) const {
	return edge_ins_cost_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Constant<UserNodeLabel, UserEdgeLabel>::
edge_del_cost_fun(const UserEdgeLabel & edge_label) const {
	return edge_del_cost_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Constant<UserNodeLabel, UserEdgeLabel>::
edge_rel_cost_fun(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const {
	if (edge_label_1 != edge_label_2) {
		return edge_rel_cost_;
	}
	return 0.0;
}

}

#endif /* SRC_EDIT_COSTS_CONSTANT_IPP_ */
