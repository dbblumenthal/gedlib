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
 * @file grec_1.ipp
 * @brief
 */

#ifndef SRC_EDIT_COSTS_GREC_1_IPP_
#define SRC_EDIT_COSTS_GREC_1_IPP_


namespace ged {

template<>
GREC1<GXLLabel, GXLLabel>::
~GREC1() {}

template<>
GREC1<GXLLabel, GXLLabel>::
GREC1() {}

template<>
double
GREC1<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return 45;
}

template<>
double
GREC1<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return 45;
}

template<>
double
GREC1<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	if (node_label_1.at("type") != node_label_2.at("type")) {
		return  90;
	}
	double x_l_minus_x_r(std::stod(node_label_1.at("x")) - std::stod(node_label_2.at("x")));
	double y_l_minus_y_r(std::stod(node_label_1.at("y")) - std::stod(node_label_2.at("y")));
	return 0.5 * std::sqrt(std::pow(x_l_minus_x_r, 2) + std::pow(y_l_minus_y_r, 2));
}

template<>
double
GREC1<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return 7.5 * std::stod(edge_label.at("frequency"));
}

template<>
double
GREC1<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return 7.5 * std::stod(edge_label.at("frequency"));
}

template<>
double
GREC1<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	if (edge_label_1.at("frequency") != edge_label_2.at("frequency")) {
		return 7.5;
	}
	if (edge_label_1.at("frequency") == "1") {
		if (edge_label_1.at("type0") == edge_label_2.at("type0")) {
			return 0;
		}
		return 15;
	}
	return 0;
}

}


#endif /* SRC_EDIT_COSTS_GREC_1_IPP_ */
