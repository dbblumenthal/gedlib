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
 * @file  ibd.ipp
 * @brief ged::IBD class definition.
 */

#ifndef SRC_EDIT_COSTS_IBD_IPP_
#define SRC_EDIT_COSTS_IBD_IPP_

namespace ged {

template<>
IBD<GXLLabel, GXLLabel>::
~IBD() {}

template<>
IBD<GXLLabel, GXLLabel>::
IBD(const std::string & otu_distances, double alpha) :
node_rel_costs_(),
alpha_{alpha},
otus_() {
	std::ifstream csv_file(otu_distances);
	std::string row;
	std::getline(csv_file, row);
	std::vector<std::string> row_as_vector;
	util::tokenize(row, ',', row_as_vector);
	std::size_t max_otu{0};
	std::vector<std::string> otu_as_vector;
	for (const std::string & otu_with_prefix : row_as_vector) {
		otu_as_vector.clear();
		util::tokenize(otu_with_prefix, '_', otu_as_vector);
		otus_.emplace_back(std::stoul(otu_as_vector.at(1)));
		max_otu = std::max(max_otu, otus_.back());
	}
	node_rel_costs_.resize(max_otu + 1, max_otu + 1);
	node_rel_costs_.set_to_val(0);

	std::size_t otu_1;
	std::size_t otu_2;
	while(std::getline(csv_file, row)) {
		row_as_vector.clear();
		util::tokenize(row, ',', row_as_vector);
		otu_as_vector.clear();
		util::tokenize(row_as_vector.at(0), '_', otu_as_vector);
		otu_1 = std::stoul(otu_as_vector.at(1));
		for (std::size_t pos{1}; pos < row_as_vector.size(); pos++) {
			otu_2 = otus_.at(pos - 1);
			node_rel_costs_(otu_1, otu_2) = std::stod(row_as_vector.at(pos));
		}
	}
	node_rel_costs_ /= node_rel_costs_.max();
}

template<>
double
IBD<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return alpha_;
}

template<>
double
IBD<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return alpha_;
}

template<>
double
IBD<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	return alpha_ * node_rel_costs_(std::stoul(node_label_1.at("OTU")), std::stoul(node_label_2.at("OTU")));
}

template<>
GXLLabel
IBD<GXLLabel, GXLLabel>::
median_node_label(const std::vector<GXLLabel> & node_labels) const {
	// Transform labels to indices.
	std::vector<std::size_t> otus;
	for (const auto & label : node_labels) {
		otus.emplace_back(std::stoul(label.at("OTU")));
	}

	// Determine the feature with the smallest SOD to all node features.
	std::size_t best_otu{undefined()};
	double best_sod{std::numeric_limits<double>::infinity()};
	double current_sod{0.0};
	for (auto otu_1 : otus_) {
		current_sod = 0.0;
		for (auto otu_2 : otus) {
			current_sod += node_rel_costs_(otu_1, otu_2);
		}
		if (current_sod < best_sod) {
			best_sod = current_sod;
			best_otu = otu_1;
		}
	}

	// Construct and return the median label.
	ged::GXLLabel median_label;
	median_label["OTU"] = std::to_string(best_otu);
	return median_label;
}

template<>
double
IBD<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_);
}

template<>
double
IBD<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_);
}

template<>
double
IBD<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	double nlogratio_1{std::stod(edge_label_1.at("nlogratio"))};
	double nlogratio_2{std::stod(edge_label_2.at("nlogratio"))};
	return (1 - alpha_) * std::fabs(nlogratio_1 - nlogratio_2);
}

template<>
GXLLabel
IBD<GXLLabel, GXLLabel>::
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

}

#endif /* SRC_EDIT_COSTS_IBD_IPP_ */
