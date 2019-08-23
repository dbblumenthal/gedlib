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
 * @file edit_costs.ipp
 * @brief ged::EditCosts class definition.
 */

#ifndef SRC_EDIT_COSTS_EDIT_COSTS_IPP_
#define SRC_EDIT_COSTS_EDIT_COSTS_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
EditCosts<UserNodeLabel, UserEdgeLabel>::
~EditCosts() {}

template<class UserNodeLabel, class UserEdgeLabel>
EditCosts<UserNodeLabel, UserEdgeLabel>::
EditCosts() {}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
node_del_cost_fun(const UserNodeLabel & node_label) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
node_ins_cost_fun(const UserNodeLabel & node_label) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
EditCosts<UserNodeLabel, UserEdgeLabel>::
vectorize_node_label(const UserNodeLabel & node_label, std::vector<double> & vector_representation) const {
	vector_representation.clear();
}

template<class UserNodeLabel, class UserEdgeLabel>
UserNodeLabel
EditCosts<UserNodeLabel, UserEdgeLabel>::
median_node_label(const std::vector<UserNodeLabel> & node_labels) const {
	return node_labels.at(0);
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
edge_rel_cost_fun(const UserEdgeLabel & node_label_1, const UserEdgeLabel & node_label_2) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
edge_del_cost_fun(const UserEdgeLabel & node_label) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
EditCosts<UserNodeLabel, UserEdgeLabel>::
edge_ins_cost_fun(const UserEdgeLabel & node_label) const {
	return 0.0;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
EditCosts<UserNodeLabel, UserEdgeLabel>::
vectorize_edge_label(const UserEdgeLabel & edge_label, std::vector<double> & vector_representation) const {
	vector_representation.clear();
}

template<class UserNodeLabel, class UserEdgeLabel>
UserEdgeLabel
EditCosts<UserNodeLabel, UserEdgeLabel>::
median_edge_label(const std::vector<UserEdgeLabel> & edge_labels) const {
	return edge_labels.at(0);
}

}

#endif /* SRC_EDIT_COSTS_EDIT_COSTS_IPP_ */
