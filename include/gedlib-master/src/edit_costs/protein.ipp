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
 * @file protein.ipp
 * @brief ged::Protein class definition.
 */

#ifndef SRC_EDIT_COSTS_PROTEIN_IPP_
#define SRC_EDIT_COSTS_PROTEIN_IPP_

namespace ged {

template<>
Protein<GXLLabel, GXLLabel>::
~Protein() {}

template<>
Protein<GXLLabel, GXLLabel>::
Protein(double node_ins_del_cost, double edge_ins_del_cost, double alpha) :
node_ins_del_cost_{node_ins_del_cost},
edge_ins_del_cost_{edge_ins_del_cost},
alpha_{alpha} {}

template<>
double
Protein<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
Protein<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * node_ins_del_cost_;
}

template<>
double
Protein<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	if (node_label_1.at("type") != node_label_2.at("type")) {
		return  alpha_ * 2 * node_ins_del_cost_;
	}
	return (alpha_ * levenshtein_distance_(node_label_1.at("sequence"), node_label_2.at("sequence")));
}

template<>
double
Protein<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * edge_ins_del_cost_ * std::stod(edge_label.at("frequency"));
}

template<>
double
Protein<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * edge_ins_del_cost_ * std::stod(edge_label.at("frequency"));
}

template<>
double
Protein<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	DMatrix edge_rel_cost_matrix(std::stod(edge_label_1.at("frequency")) + 1, std::stod(edge_label_2.at("frequency")) + 1);
	for (std::size_t i{0}; i < edge_rel_cost_matrix.num_rows() - 1; i++) {
		edge_rel_cost_matrix(i, edge_rel_cost_matrix.num_cols() - 1) = edge_ins_del_cost_;
	}
	for (std::size_t k{0}; k < edge_rel_cost_matrix.num_cols() - 1; k++) {
		edge_rel_cost_matrix(edge_rel_cost_matrix.num_rows() - 1, k) = edge_ins_del_cost_;
	}
	for (std::size_t i{0}; i < edge_rel_cost_matrix.num_rows() - 1; i++) {
		for (std::size_t k{0}; k < edge_rel_cost_matrix.num_cols() - 1; k++) {
			if (edge_label_1.at(std::string("type") + std::to_string(i)) == edge_label_2.at(std::string("type") + std::to_string(k))) {
				edge_rel_cost_matrix(i, k) = 0;
			}
			else {
				edge_rel_cost_matrix(i, k) = 2 * edge_ins_del_cost_;
			}
		}
	}
	LSAPESolver lsape_solver(&edge_rel_cost_matrix);
	lsape_solver.solve();
	return (1 - alpha_) * lsape_solver.minimal_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Protein<UserNodeLabel, UserEdgeLabel>::
levenshtein_distance_(const std::string & string_1, const std::string & string_2) const {
	DMatrix prefix_distance(string_1.size() + 1, string_2.size() + 1);
	prefix_distance(0, 0) = 0.0;
	for (std::size_t prefix_size_1{1}; prefix_size_1 <= string_1.size(); prefix_size_1++) {
		prefix_distance(prefix_size_1, 0) = static_cast<double>(prefix_size_1);
	}
	for (std::size_t prefix_size_2{1}; prefix_size_2 <= string_2.size(); prefix_size_2++) {
		prefix_distance(0, prefix_size_2) = static_cast<double>(prefix_size_2);
	}
	for (std::size_t prefix_size_1{1}; prefix_size_1 <= string_1.size(); prefix_size_1++) {
		for (std::size_t prefix_size_2{1}; prefix_size_2 <= string_2.size(); prefix_size_2++) {
			prefix_distance(prefix_size_1, prefix_size_2) = prefix_distance(prefix_size_1 - 1, prefix_size_2 - 1);
			if (string_1.at(prefix_size_1 - 1) != string_2.at(prefix_size_2 - 1)) {
				prefix_distance(prefix_size_1, prefix_size_2) += 1;
			}
			prefix_distance(prefix_size_1, prefix_size_2) = std::min(prefix_distance(prefix_size_1, prefix_size_2), prefix_distance(prefix_size_1 - 1, prefix_size_2) + 1);
			prefix_distance(prefix_size_1, prefix_size_2) = std::min(prefix_distance(prefix_size_1, prefix_size_2), prefix_distance(prefix_size_1, prefix_size_2 - 1) + 1);
		}
	}
	return prefix_distance(string_1.size(), string_2.size());
}

}

#endif /* SRC_EDIT_COSTS_PROTEIN_IPP_ */
