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
 * @file chem_1.ipp
 * @brief ged::CHEM1 class definition.
 */

#ifndef SRC_EDIT_COSTS_CHEM_1_IPP_
#define SRC_EDIT_COSTS_CHEM_1_IPP_

namespace ged {

template<>
CHEM1<GXLLabel, GXLLabel>::
~CHEM1() {}

template<>
CHEM1<GXLLabel, GXLLabel>::
CHEM1(double node_ins_del_cost, double node_rel_cost, double edge_ins_del_cost, double edge_rel_cost) :
node_ins_del_cost_{node_ins_del_cost},
node_rel_cost_{node_rel_cost},
edge_ins_del_cost_{edge_ins_del_cost},
edge_rel_cost_{edge_rel_cost} {}

template<>
double
CHEM1<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return node_ins_del_cost_;
}

template<>
double
CHEM1<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return node_ins_del_cost_;
}

template<>
double
CHEM1<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	if (node_label_1.at("chem") != node_label_2.at("chem")) {
		return node_rel_cost_;
	}
	return 0.0;
}

template<>
double
CHEM1<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return edge_ins_del_cost_;
}

template<>
double
CHEM1<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return edge_ins_del_cost_;
}

template<>
double
CHEM1<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	if (edge_label_1.at("valence") != edge_label_2.at("valence")) {
		return edge_rel_cost_;
	}
	return 0.0;
}

}

#endif /* SRC_EDIT_COSTS_CHEM_1_IPP_ */
