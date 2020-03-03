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
IBDCosts(const std::string & distance_matrix, double alpha, const std::string & feature_name) :
node_rel_costs_(0, 0),
feature_name_(feature_name),
alpha_{alpha},
feature_id_to_index_(),
index_to_feature_id_() {
	std::ifstream csv_file(distance_matrix);
	std::string row;
	std::getline(csv_file, row);
	std::vector<std::string> row_as_vector;
	util::tokenize(row, ',', row_as_vector);
	std::size_t num_features{row_as_vector.size() - 1};
	node_rel_costs_.resize(num_features, num_features);
	std::size_t feature_id;
	std::size_t max_feature_id{0};
	for (std::size_t pos{1}; pos < row_as_vector.size(); pos++) {
		feature_id = std::stoul(row_as_vector.at(pos));
		index_to_feature_id_.emplace_back(feature_id);
		max_feature_id = std::max(max_feature_id, feature_id);
	}
	feature_id_to_index_ = std::vector<std::size_t>(max_feature_id + 1, undefined());
	for (std::size_t index{0}; index < num_features; index++) {
		feature_id_to_index_.at(index_to_feature_id_.at(index)) = index;
	}
	std::size_t index_1{0};
	while(std::getline(csv_file, row)) {
		row_as_vector.clear();
		util::tokenize(row, ',', row_as_vector);
		feature_id = std::stoul(row_as_vector.at(0));
		for (std::size_t index_2{0}; index_2 < num_features; index_2++) {
			node_rel_costs_(index_1, index_2) = std::stod(row_as_vector.at(index_2 + 1));
		}
		index_1++;
	}
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return alpha_;
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return alpha_;
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	std::size_t feature_id_1{std::stoul(node_label_1.at(feature_name_))};
	std::size_t feature_id_2{std::stoul(node_label_2.at(feature_name_))};
	return alpha_ * node_rel_costs_(feature_id_to_index_.at(feature_id_1), feature_id_to_index_.at(feature_id_2));
}

template<>
GXLLabel
IBDCosts<GXLLabel, GXLLabel>::
median_node_label(const std::vector<GXLLabel> & node_labels) const {
	// Transform labels to indices.
	std::vector<std::size_t> indices;
	for (const auto & label : node_labels) {
		indices.emplace_back(feature_id_to_index_.at(std::stoul(label.at(feature_name_))));
	}

	// Determine the feature with the smallest SOD to all node features.
	std::size_t best_index{undefined()};
	double best_sod{std::numeric_limits<double>::infinity()};
	double current_sod{0.0};
	for (std::size_t index_1{0}; index_1 < index_to_feature_id_.size(); index_1++) {
		current_sod = 0.0;
		for (std::size_t index_2 : indices) {
			current_sod += node_rel_costs_(index_1, index_2);
		}
		if (current_sod < best_sod) {
			best_sod = current_sod;
			best_index = index_1;
		}
	}

	// Construct and return the median label.
	ged::GXLLabel median_label;
	median_label[feature_name_] = std::to_string(index_to_feature_id_.at(best_index));
	return median_label;
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_);
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	return (1 - alpha_);
}

template<>
double
IBDCosts<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	double nlogratio_1{std::stod(edge_label_1.at("nlogratio"))};
	double nlogratio_2{std::stod(edge_label_2.at("nlogratio"))};
	return (1 - alpha_) * std::fabs(nlogratio_1 - nlogratio_2);
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

#endif /* APP_IBD_SRC_IBD_COSTS_IPP_ */
