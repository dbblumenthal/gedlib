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
 * @file cmu.ipp
 * @brief ged::CMU class definition.
 */

#ifndef SRC_EDIT_COSTS_CMU_IPP_
#define SRC_EDIT_COSTS_CMU_IPP_

namespace ged {

template<>
CMU<GXLLabel, GXLLabel>::
~CMU() {}

template<>
CMU<GXLLabel, GXLLabel>::
CMU(double node_ins_del_cost, double alpha) :
node_ins_del_cost_{node_ins_del_cost},
alpha_{alpha} {}

template<>
double
CMU<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
CMU<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
CMU<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	return 0.0;
}

template<>
void
CMU<GXLLabel, GXLLabel>::
vectorize_node_label(const GXLLabel & node_label, std::vector<double> & vector_representation) const {
	vector_representation.clear();
	vector_representation.push_back(std::stod(node_label.at("x")));
	vector_representation.push_back(std::stod(node_label.at("y")));
}

template<>
double
CMU<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * std::stod(edge_label.at("dist"));
}

template<>
double
CMU<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * std::stod(edge_label.at("dist"));
}

template<>
double
CMU<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	return (1 - alpha_) * std::abs(std::stod(edge_label_1.at("dist")) - std::stod(edge_label_2.at("dist")));
}

}

#endif /* SRC_EDIT_COSTS_CMU_IPP_ */
