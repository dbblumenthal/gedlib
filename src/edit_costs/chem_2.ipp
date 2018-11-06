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
 * @file chem_2.ipp
 * @brief ged::CHEM2 class definition.
 */

#ifndef SRC_EDIT_COSTS_CHEM_2_IPP_
#define SRC_EDIT_COSTS_CHEM_2_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
CHEM2<UserNodeLabel, UserEdgeLabel>::
~CHEM2() {}

template<class UserNodeLabel, class UserEdgeLabel>
CHEM2<UserNodeLabel, UserEdgeLabel>::
CHEM2(double node_ins_del_cost, double edge_ins_del_cost, double alpha) :
node_ins_del_cost_{node_ins_del_cost},
edge_ins_del_cost_{edge_ins_del_cost},
alpha_{alpha} {}

template<class UserNodeLabel, class UserEdgeLabel>
double
CHEM2<UserNodeLabel, UserEdgeLabel>::
node_ins_cost_fun(const UserNodeLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
CHEM2<UserNodeLabel, UserEdgeLabel>::
node_del_cost_fun(const UserNodeLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	if (node_label_1.at("chem") != node_label_2.at("chem")) {
		return alpha_ * 2 * node_ins_del_cost_;
	}
	return 0.0;
}

template<>
double
CHEM2<GXLLabel, double>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	if (node_label_1.at("chem") != node_label_2.at("chem")) {
		return alpha_ * 2 * node_ins_del_cost_;
	}
	return 0.0;
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * edge_ins_del_cost_;
}

template<>
double
CHEM2<GXLLabel, double>::
edge_ins_cost_fun(const double & edge_label) const {
	return (1 - alpha_) * std::fabs(edge_label);
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * edge_ins_del_cost_;
}

template<>
double
CHEM2<GXLLabel, double>::
edge_del_cost_fun(const double & edge_label) const {
	return (1 - alpha_) * std::fabs(edge_label);
}

template<>
double
CHEM2<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	if (edge_label_1.at("valence") != edge_label_2.at("valence")) {
		return (1 - alpha_) * edge_ins_del_cost_;
	}
	return 0.0;
}

template<>
double
CHEM2<GXLLabel, double>::
edge_rel_cost_fun(const double & edge_label_1, const double & edge_label_2) const {
	return (1 - alpha_) * std::fabs(edge_label_1 - edge_label_2);
}

}




#endif /* SRC_EDIT_COSTS_CHEM_2_IPP_ */
