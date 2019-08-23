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
 * @file ring_ml.ipp
 * @brief RingML class definition.
 */

#ifndef SRC_METHODS_RING_ML_IPP_
#define SRC_METHODS_RING_ML_IPP_

namespace ged {

// === Destructors and constructors. ===
template<class UserNodeLabel, class UserEdgeLabel>
RingML<UserNodeLabel, UserEdgeLabel>::
~RingML() {}

template<class UserNodeLabel, class UserEdgeLabel>
RingML<UserNodeLabel, UserEdgeLabel>::
RingML(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
MLBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
rings_(),
led_method_{LSAPE_OPTIMAL},
sort_method_{COUNTING},
use_topological_features_{true},
use_global_features_{true},
num_layers_{undefined()},
global_features_() {}

// === Definitions of member functions inherited from MLBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
ml_init_graph_(const GEDGraph & graph) {
	build_rings_(graph);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
ml_set_default_options_() {
	led_method_ = LSAPE_OPTIMAL;
	sort_method_ = COUNTING;
	use_topological_features_ = true;
	use_global_features_ = true;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
RingML<UserNodeLabel, UserEdgeLabel>::
ml_valid_options_string_() const {
	return "[--led-method <arg>] [--sort-method <arg>] [--topological-features <arg>] [--global-features <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
RingML<UserNodeLabel, UserEdgeLabel>::
ml_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "led-method") {
		if (arg == "LSAPE_OPTIMAL") {
			led_method_ = LSAPE_OPTIMAL;
		}
		else if (arg  == "LSAPE_GREEDY") {
			led_method_ = LSAPE_GREEDY;
		}
		else if (arg == "GAMMA") {
			led_method_ = GAMMA;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option led-method. Usage: options = \"[--led-method LSAPE_OPTIMAL|LSAPE_GREEDY|GAMMA] [...]\"");
		}
		return true;
	}
	else if (option == "sort-method") {
		if (arg == "COUNTING") {
			sort_method_ = COUNTING;
		}
		else if (arg  == "STD") {
			sort_method_ = STD;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option sort-method. Usage: options = \"[--sort-method COUNTING|STD] [...]\"");
		}
		return true;
	}
	else if (option == "topological-features") {
		if (arg == "TRUE") {
			use_topological_features_ = true;
		}
		else if (arg  == "FALSE") {
			use_topological_features_ = false;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option topological-features. Usage: options = \"[--topological-features TRUE|FALSE] [...]\"");
		}
		return true;
	}
	else if (option == "global-features") {
		if (arg == "TRUE") {
			use_global_features_ = true;
		}
		else if (arg  == "FALSE") {
			use_global_features_ = false;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option global-features. Usage: options = \"[--global-features TRUE|FALSE] [...]\"");
		}
		return true;
	}
	return false;

}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
ml_init_feature_variables_(const GEDGraph & g, const GEDGraph & h, std::size_t num_threads) {
	if (use_global_features_) {
		global_features_.clear();
		global_features_.push_back(static_cast<double>(g.num_nodes()));
		global_features_.push_back(static_cast<double>(h.num_nodes()));
		global_features_.push_back(this->ged_data_.mean_node_subs_cost(g, h));
		global_features_.push_back(this->ged_data_.mean_node_del_cost(g));
		global_features_.push_back(this->ged_data_.mean_node_ins_cost(h));
		global_features_.push_back(static_cast<double>(g.num_edges()));
		global_features_.push_back(static_cast<double>(h.num_edges()));
		global_features_.push_back(this->ged_data_.mean_edge_subs_cost(g, h));
		global_features_.push_back(this->ged_data_.mean_edge_del_cost(g));
		global_features_.push_back(this->ged_data_.mean_edge_ins_cost(h));
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
ml_populate_substitution_feature_vector_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k, std::vector<double> & feature_vector) {
	feature_vector.clear();
	add_global_features_(feature_vector);
	const Ring_ & ring_i = rings_.at(g.id()).at(i);
	const Ring_ & ring_k = rings_.at(h.id()).at(k);
	for (std::size_t level{0}; level < num_layers_; level++) {
		add_layer_substitution_features_(ring_i, ring_k, level, feature_vector);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
ml_populate_deletion_feature_vector_(const GEDGraph & g, GEDGraph::NodeID i, std::vector<double> & feature_vector) {
	feature_vector.clear();
	add_global_features_(feature_vector);
	const Ring_ & ring_i = rings_.at(g.id()).at(i);
	for (std::size_t level{0}; level < num_layers_; level++) {
		add_layer_deletion_features_(ring_i, level, feature_vector);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
ml_populate_insertion_feature_vector_(const GEDGraph & h, GEDGraph::NodeID k, std::vector<double> & feature_vector) {
	feature_vector.clear();
	add_global_features_(feature_vector);
	const Ring_ & ring_k = rings_.at(h.id()).at(k);
	for (std::size_t level{0}; level < num_layers_; level++) {
		add_layer_deletion_features_(ring_k, level, feature_vector);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
ml_init_for_num_features_() {
	if (this->num_features_ == undefined()) {
		return;
	}
	std::size_t num_local_features{this->num_features_};
	if (use_global_features_) {
		num_local_features -= 10;
	}
	if (use_topological_features_) {
		num_layers_ = num_local_features / 6;
	}
	else {
		num_layers_ = num_local_features / 3;
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
RingML<UserNodeLabel, UserEdgeLabel>::
ml_get_num_features_() {
	set_num_layers_();
	std::size_t num_features{num_layers_};
	if (use_topological_features_) {
		num_features *= 6;
	}
	else {
		num_features *= 3;
	}
	if (use_global_features_) {
		num_features += 10;
	}
	return num_features;
}

// === Defintion of private struct Layer_. ===
template<class UserNodeLabel, class UserEdgeLabel>
RingML<UserNodeLabel, UserEdgeLabel>::
Layer_ ::
Layer_(std::size_t level) :
level{level},
node_labels(),
inner_edge_labels(),
outer_edge_labels() {}

// === Definition of private struct Ring_. ===
template<class UserNodeLabel, class UserEdgeLabel>
RingML<UserNodeLabel, UserEdgeLabel>::
Ring_ ::
Ring_() :
layers() {}

// === Definitions of helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
build_rings_(const GEDGraph & graph) {
	rings_[graph.id()] = NodeRingMap_();
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		build_ring_(graph, *node, rings_.at(graph.id()));
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
build_ring_(const GEDGraph & graph, GEDGraph::NodeID root, NodeRingMap_ & rings) {
	std::map<GEDGraph::NodeID, int> distance_to_root;
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		distance_to_root[*node] = -1;
	}
	distance_to_root[root] = 0;

	std::map<GEDGraph::EdgeID, bool> discovered_edge;
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		discovered_edge[*edge] = false;
	}

	Layer_ current_layer(0);
	std::queue<GEDGraph::NodeID> queue;
	queue.push(root);
	while (not queue.empty()) {
		GEDGraph::NodeID current_node{queue.front()};
		queue.pop();
		if (static_cast<int>(current_layer.level) < distance_to_root.at(current_node)) {
			if (led_method_ == GAMMA) {
				if (sort_method_ == COUNTING) {
					util::counting_sort(current_layer.node_labels.begin(), current_layer.node_labels.end());
					util::counting_sort(current_layer.inner_edge_labels.begin(), current_layer.inner_edge_labels.end());
					util::counting_sort(current_layer.outer_edge_labels.begin(), current_layer.outer_edge_labels.end());
				}
				else {
					std::sort(current_layer.node_labels.begin(), current_layer.node_labels.end());
					std::sort(current_layer.inner_edge_labels.begin(), current_layer.inner_edge_labels.end());
					std::sort(current_layer.outer_edge_labels.begin(), current_layer.outer_edge_labels.end());
				}
			}
			rings[root].layers.push_back(current_layer);
			current_layer = Layer_(current_layer.level + 1);
		}
		current_layer.node_labels.push_back(graph.get_node_label(current_node));
		for (auto edge = graph.incident_edges(current_node).first; edge != graph.incident_edges(current_node).second; edge++) {
			GEDGraph::NodeID next_node{graph.head(*edge)};
			if (distance_to_root.at(next_node) == -1) {
				distance_to_root[next_node] = current_layer.level + 1;
				if (current_layer.level < num_layers_) {
					queue.push(next_node);
				}
			}
			if (not discovered_edge.at(*edge)) {
				discovered_edge[*edge] = true;
				if (distance_to_root.at(current_node) == distance_to_root.at(next_node)) {
					current_layer.inner_edge_labels.push_back(graph.get_edge_label(*edge));
				}
				else if (distance_to_root.at(current_node) < distance_to_root.at(next_node)) {
					current_layer.outer_edge_labels.push_back(graph.get_edge_label(*edge));
				}
				else {
					throw Error(std::string("Error when building ring rooted at ") + std::to_string(root) +
							" for graph " + std::to_string(graph.id()) + ": dist(" +
							std::to_string(current_node) +") = " + std::to_string(distance_to_root.at(current_node)) +
							" > dist(" + std::to_string(next_node) +") = " + std::to_string(distance_to_root.at(next_node)));
				}
			}
		}
	}
	if (led_method_ == GAMMA) {
		if (sort_method_ == COUNTING) {
			util::counting_sort(current_layer.node_labels.begin(), current_layer.node_labels.end());
			util::counting_sort(current_layer.inner_edge_labels.begin(), current_layer.inner_edge_labels.end());
			util::counting_sort(current_layer.outer_edge_labels.begin(), current_layer.outer_edge_labels.end());
		}
		else {
			std::sort(current_layer.node_labels.begin(), current_layer.node_labels.end());
			std::sort(current_layer.inner_edge_labels.begin(), current_layer.inner_edge_labels.end());
			std::sort(current_layer.outer_edge_labels.begin(), current_layer.outer_edge_labels.end());
		}
	}
	rings[root].layers.push_back(current_layer);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
set_num_layers_() {
	num_layers_ = 0;
	for (auto graph_rings_pair = rings_.begin(); graph_rings_pair != rings_.end(); graph_rings_pair++) {
		for (auto ring = graph_rings_pair->second.begin(); ring != graph_rings_pair->second.end(); ring++) {
			num_layers_ = std::max(num_layers_, ring->second.layers.size());
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
add_global_features_(std::vector<double> & feature_vector) const {
	for (auto feature = global_features_.begin(); feature != global_features_.end(); feature++) {
		feature_vector.push_back(*feature);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
add_layer_substitution_features_(const Ring_ & ring_i, const Ring_ & ring_k, std::size_t level, std::vector<double> & feature_vector) const {
	if ((ring_i.layers.size() > level) and (ring_k.layers.size() > level)) {
		add_layer_features_(ring_i.layers.at(level), ring_k.layers.at(level), feature_vector);
	}
	else if (ring_i.layers.size() > level) {
		add_layer_features_(ring_i.layers.at(level), Layer_(0), feature_vector);
	}
	else if (ring_k.layers.size() > level) {
		add_layer_features_(Layer_(0), ring_k.layers.at(level), feature_vector);
	}
	else {
		std::size_t num_layer_features = use_topological_features_ ? 6 : 3;
		for (std::size_t layer_feature{0}; layer_feature < num_layer_features; layer_feature++) {
			feature_vector.push_back(0.0);
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
add_layer_deletion_features_(const Ring_ & ring, std::size_t level, std::vector<double> & feature_vector) const {
	if (ring.layers.size() > level) {
		add_layer_features_(ring.layers.at(level), Layer_(0), feature_vector);
	}
	else {
		std::size_t num_layer_features = use_topological_features_ ? 6 : 3;
		for (std::size_t layer_feature{0}; layer_feature < num_layer_features; layer_feature++) {
			feature_vector.push_back(0.0);
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
add_layer_insertion_features_(const Ring_ & ring, std::size_t level, std::vector<double> & feature_vector) const {
	if (ring.layers.size() > level) {
		add_layer_features_(Layer_(0), ring.layers.at(level), feature_vector);
	}
	else {
		std::size_t num_layer_features = use_topological_features_ ? 6 : 3;
		for (std::size_t layer_feature{0}; layer_feature < num_layer_features; layer_feature++) {
			feature_vector.push_back(0.0);
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
RingML<UserNodeLabel, UserEdgeLabel>::
add_layer_features_(const Layer_ & lhs, const Layer_ & rhs, std::vector<double> & feature_vector) const {

	if (use_topological_features_) {
		feature_vector.push_back(static_cast<double>(lhs.node_labels.size() - rhs.node_labels.size()));
		feature_vector.push_back(static_cast<double>(lhs.inner_edge_labels.size() - rhs.inner_edge_labels.size()));
		feature_vector.push_back(static_cast<double>(lhs.outer_edge_labels.size() - rhs.outer_edge_labels.size()));
	}

	switch (led_method_) {
	case GAMMA:
		feature_vector.push_back(gamma_multiset_cost_(lhs.node_labels, rhs.node_labels, true));
		feature_vector.push_back(gamma_multiset_cost_(lhs.inner_edge_labels, rhs.inner_edge_labels, false));
		feature_vector.push_back(gamma_multiset_cost_(lhs.outer_edge_labels, rhs.outer_edge_labels, false));
		break;
	default:
		feature_vector.push_back(lsape_multiset_cost_(lhs.node_labels, rhs.node_labels, true));
		feature_vector.push_back(lsape_multiset_cost_(lhs.inner_edge_labels, rhs.inner_edge_labels, false));
		feature_vector.push_back(lsape_multiset_cost_(lhs.outer_edge_labels, rhs.outer_edge_labels, false));
		break;
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
RingML<UserNodeLabel, UserEdgeLabel>::
lsape_multiset_cost_(const std::vector<LabelID> & lhs, const std::vector<LabelID> & rhs, bool node_labels) const {

	if ((lhs.size() == 0) and (rhs.size() == 0)) {
		return 0.0;
	}

	if ((lhs.size() > 0) and (rhs.size() == 0)) {
		double cost{0.0};
		for (std::size_t row{0}; row < lhs.size(); row++) {
			if (node_labels) {
				cost += this->ged_data_.node_cost(lhs.at(row), dummy_label());
			}
			else {
				cost += this->ged_data_.edge_cost(lhs.at(row), dummy_label());
			}
		}
		return cost;
	}

	if ((lhs.size() == 0) and (rhs.size() > 0)) {
		double cost{0.0};
		for (std::size_t col{0}; col < rhs.size(); col++) {
			if (node_labels) {
				cost += this->ged_data_.node_cost(dummy_label(), rhs.at(col));
			}
			else {
				cost += this->ged_data_.edge_cost(dummy_label(), rhs.at(col));
			}
		}

		return cost;
	}

	DMatrix problem(lhs.size() + 1, rhs.size() + 1, 0.0);
	// Collect deletion costs.
	for (std::size_t row{0}; row < lhs.size(); row++) {
		if (node_labels) {
			problem(row, rhs.size()) = this->ged_data_.node_cost(lhs.at(row), dummy_label());
		}
		else {
			problem(row, rhs.size()) = this->ged_data_.edge_cost(lhs.at(row), dummy_label());
		}
	}

	// Collect insertion costs.
	for (std::size_t col{0}; col < rhs.size(); col++) {
		if (node_labels) {
			problem(lhs.size(), col) = this->ged_data_.node_cost(dummy_label(), rhs.at(col));
		}
		else {
			problem(lhs.size(), col) = this->ged_data_.edge_cost(dummy_label(), rhs.at(col));
		}
	}

	// Collect substitution costs.
	for (std::size_t row{0}; row < lhs.size(); row++) {
		for (std::size_t col{0}; col < rhs.size(); col++) {
			if (node_labels) {
				problem(row, col) = this->ged_data_.node_cost(lhs.at(row), rhs.at(col));
			}
			else {
				problem(row, col) = this->ged_data_.edge_cost(lhs.at(row), rhs.at(col));
			}
		}
	}

	LSAPESolver problem_solver(&problem);
	if (led_method_ == LSAPE_OPTIMAL) {
		problem_solver.set_model(this->lsape_model_);
	}
	else {
		problem_solver.set_greedy_method(this->greedy_method_);
	}
	problem_solver.solve();

	return problem_solver.minimal_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
double
RingML<UserNodeLabel, UserEdgeLabel>::
gamma_multiset_cost_(const std::vector<LabelID> & lhs, const std::vector<LabelID> & rhs, bool node_labels) const {
	double cost{0.0};

	// Compute and add minimal edge insertion costs.
	if (lhs.size() < rhs.size()) {
		double avg_ins_cost{0.0};
		for (auto label_rhs = rhs.begin(); label_rhs != rhs.end(); label_rhs++) {
			if (node_labels) {
				avg_ins_cost += this->ged_data_.node_cost(dummy_label(), *label_rhs);
			}
			else {
				avg_ins_cost += this->ged_data_.edge_cost(dummy_label(), *label_rhs);
			}
		}
		avg_ins_cost /= static_cast<double>(rhs.size());
		cost += static_cast<double>(rhs.size() - lhs.size()) * avg_ins_cost;
	}

	// Compute and add minimal edge deletion costs.
	if (lhs.size() > rhs.size()) {
		double avg_del_cost{0.0};
		for (auto label_lhs = lhs.begin(); label_lhs != lhs.end(); label_lhs++) {
			if (node_labels) {
				avg_del_cost += this->ged_data_.node_cost(*label_lhs, dummy_label());
			}
			else {
				avg_del_cost += this->ged_data_.edge_cost(*label_lhs, dummy_label());
			}
		}
		avg_del_cost /= static_cast<double>(lhs.size());
		cost += static_cast<double>(lhs.size() - rhs.size()) * avg_del_cost;
	}

	// Compute minimal edge relabelling costs.
	double avg_rel_cost{0.0};
	std::size_t count_diff_labels{0};
	for (auto label_lhs = lhs.begin(); label_lhs != lhs.end(); label_lhs++) {
		for (auto label_rhs = rhs.begin(); label_rhs != rhs.end(); label_rhs++) {
			if (*label_lhs != *label_rhs) {
				count_diff_labels++;
				if (node_labels) {
					avg_rel_cost += this->ged_data_.node_cost(*label_lhs, *label_rhs);
				}
				else {
					avg_rel_cost += this->ged_data_.edge_cost(*label_lhs, *label_rhs);
				}
			}
		}
	}
	avg_rel_cost /= static_cast<double>(count_diff_labels);

	// Compute multiset intersection size.
	std::size_t intersection_size{0};
	auto label_lhs = lhs.begin();
	auto label_rhs = rhs.begin();
	while ((label_lhs != lhs.end()) and (label_rhs != rhs.end())) {
		if (*label_lhs == *label_rhs) {
			intersection_size++;
			label_lhs++;
			label_rhs++;
		}
		else if (*label_lhs < *label_rhs) {
			label_lhs++;
		}
		else {
			label_rhs++;
		}
	}

	std::size_t gamma(std::min(lhs.size(), rhs.size()) - intersection_size);
	if (gamma > 0) {
		cost += static_cast<double>(gamma) * avg_rel_cost;
	}

	return cost;
}

}

#endif /* SRC_METHODS_RING_ML_IPP_ */
