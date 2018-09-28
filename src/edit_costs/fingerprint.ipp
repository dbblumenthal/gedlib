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
 * @file fingerprint.ipp
 * @brief ged::Fingerprint class definition.
 */

#ifndef SRC_EDIT_COSTS_FINGERPRINT_IPP_
#define SRC_EDIT_COSTS_FINGERPRINT_IPP_


namespace ged {

template<>
Fingerprint<GXLLabel, GXLLabel>::
~Fingerprint() {}

template<>
Fingerprint<GXLLabel, GXLLabel>::
Fingerprint(double node_ins_del_cost, double edge_ins_del_cost, double alpha) :
node_ins_del_cost_{node_ins_del_cost},
edge_ins_del_cost_{edge_ins_del_cost},
alpha_{alpha} {}

template<>
double
Fingerprint<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
Fingerprint<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
Fingerprint<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	return 0;
}

template<>
double
Fingerprint<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * edge_ins_del_cost_;
}

template<>
double
Fingerprint<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * edge_ins_del_cost_;
}

template<>
double
Fingerprint<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	double orient_diff{std::fabs(std::stod(edge_label_1.at("orient")) - std::stod(edge_label_2.at("orient")))};
	return (1 - alpha_) * std::min(pi_() - orient_diff, orient_diff);
}

}


#endif /* SRC_EDIT_COSTS_FINGERPRINT_IPP_ */
