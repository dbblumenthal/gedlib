/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2020 by David B. Blumenthal                              *
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
 * @file  ibd_costs.ipp
 * @brief ged::IBDCosts class definition.
 */

#ifndef APP_IBD_SRC_IBD_COSTS_IPP_
#define APP_IBD_SRC_IBD_COSTS_IPP_

using namespace ged;

template<>
IBDCosts<GXLLabel, GXLLabel>::
~IBDCosts() {}

template<>
IBDCosts<GXLLabel, GXLLabel>::
IBDCosts(const std::string & distance_matrix, double ins_del_factor, double alpha) :
node_rel_costs_(452, 452),
max_node_rel_cost_{0.0},
otu_to_index_(2197, undefined()),
index_to_otu_(452, undefined()),
max_edge_rel_cost_{1.0},
alpha_{alpha},
ins_del_factor_{ins_del_factor} {
	std::ifstream csv_file(distance_matrix);
	std::string row;
	std::getline(csv_file, row);
	std::size_t max_otu_index{451};
	std::size_t max_otu{2196};
	std::vector<std::string> row_as_vector;
	std::size_t otu;
	std::size_t otu_index_1{0};
	while(std::getline(csv_file, row)) {
		row_as_vector.clear();
		util::tokenize(row, ',', row_as_vector);
		if (row_as_vector.size() != max_otu_index + 2) {
			throw Error("Each row in the distance matrix must have exactly 453 columns.");
		}
		otu = std::stoul(row_as_vector.at(0));
		if (otu > max_otu) {
			throw Error(std::string("OTU ") + std::to_string(otu) + " is larger than 2196.");
		}
		otu_to_index_.at(otu) = otu_index_1;
		index_to_otu_.at(otu_index_1) = otu;
		for (std::size_t otu_index_2{0}; otu_index_2 <= max_otu_index; otu_index_2++) {
			node_rel_costs_(otu_index_1, otu_index_2) = std::stod(row_as_vector.at(otu_index_2 + 1));
		}
		otu_index_1++;
	}
	if (otu_index_1 != max_otu_index + 1) {
		throw Error("The distance matrix must have exactly 453 rows.");
	}
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * ins_del_factor_ * max_node_rel_cost_;
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return alpha_ * ins_del_factor_ * max_node_rel_cost_;
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	std::size_t otu_1{std::stoul(node_label_1.at("OTU"))};
	std::size_t otu_2{std::stoul(node_label_2.at("OTU"))};
	return node_rel_costs_(otu_to_index_.at(otu_1), otu_to_index_.at(otu_2));
}

template<>
GXLLabel
IBDCosts<GXLLabel, GXLLabel>::
median_node_label(const std::vector<GXLLabel> & node_labels) const {
	// Transform labels to OTU indices.
	std::vector<std::size_t> otu_indices;
	for (const auto & label : node_labels) {
		otu_indices.emplace_back(otu_to_index_.at(std::stoul(label.at("OTU"))));
	}

	// Determine the OTU with the smallest SOD to the OTUs of the node labels.
	std::size_t best_otu_index{undefined()};
	double best_sod{std::numeric_limits<double>::infinity()};
	double current_sod{0.0};
	for (std::size_t otu_index_1{0}; otu_index_1 < index_to_otu_.size(); otu_index_1++) {
		current_sod = 0.0;
		for (std::size_t otu_index_2 : otu_indices) {
			current_sod += node_rel_costs_(otu_index_1, otu_index_2);
		}
		if (current_sod < best_sod) {
			best_sod = current_sod;
			best_otu_index = otu_index_1;
		}
	}

	// Construct and return the median label.
	ged::GXLLabel median_label;
	median_label["OTU"] = std::to_string(index_to_otu_.at(best_otu_index));
	return median_label;
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * ins_del_factor_ * max_edge_rel_cost_;
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_) * ins_del_factor_ * max_edge_rel_cost_;
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	double nlogratio_1{std::stod(edge_label_1.at("nlogratio"))};
	double nlogratio_2{std::stod(edge_label_2.at("nlogratio"))};
	return std::fabs(nlogratio_1 - nlogratio_2);
}

template<>
GXLLabel
IBDCosts<GXLLabel, GXLLabel>::
median_edge_label(const std::vector<GXLLabel> & edge_labels) const {

	// Collect normalized log-ratios.
	std::vector<double> nlogratios;
	for (const auto & label : edge_labels) {
		nlogratios.emplace_back(std::stod(label.at("nlogratio")));
	}

	// Compute the median.
	double median;
	const auto middle_itr = nlogratios.begin() + nlogratios.size() / 2;
	std::nth_element(nlogratios.begin(), middle_itr, nlogratios.end());
	if (nlogratios.size() % 2 == 0) {
		const auto left_middle_itr = std::max_element(nlogratios.begin(), middle_itr);
		median = (*left_middle_itr + *middle_itr) / 2.0;
	}
	else {
		median = *middle_itr;
	}

	// Construct and return the median label.
	ged::GXLLabel median_label;
	median_label["nlogratio"] = std::to_string(median);
	return median_label;
}

template<>
void
IBDCosts<GXLLabel, GXLLabel>::
set_alpha(double alpha) {
	alpha_ = alpha;
}

template<>
void
IBDCosts<GXLLabel, GXLLabel>::
set_ins_del_factor(double ins_del_factor) {
	ins_del_factor_ = ins_del_factor;
}

#endif /* APP_IBD_SRC_IBD_COSTS_IPP_ */
