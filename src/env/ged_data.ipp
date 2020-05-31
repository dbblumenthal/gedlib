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
 * @file  ged_data.ipp
 * @brief ged::GEDData class definition.
 */

#ifndef SRC_ENV_GED_DATA_IPP_
#define SRC_ENV_GED_DATA_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
GEDData<UserNodeLabel, UserEdgeLabel>::
GEDData() :
graphs_(),
graph_names_(),
graph_classes_(),
num_graphs_without_shuffled_copies_{0},
strings_to_internal_node_ids_(),
internal_node_ids_to_strings_(),
edit_costs_{nullptr},
node_costs_(),
edge_costs_(),
node_labels_(),
node_label_ids_(),
edge_labels_(),
edge_label_ids_(),
init_type_{Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES},
delete_edit_costs_{true},
max_num_nodes_{0},
max_num_edges_{0} {}

template<class UserNodeLabel, class UserEdgeLabel>
GEDData<UserNodeLabel, UserEdgeLabel>::
~GEDData() {
	if (delete_edit_costs_) {
		delete edit_costs_;
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
set_edit_costs_(Options::EditCosts edit_costs, const std::vector<double> & edit_cost_constants) {
	if (delete_edit_costs_) {
		delete edit_costs_;
	}
	switch (edit_costs) {
	case Options::EditCosts::CONSTANT:
		if (edit_cost_constants.size() == 6) {
			edit_costs_ = new Constant<UserNodeLabel, UserEdgeLabel>(edit_cost_constants.at(0), edit_cost_constants.at(1), edit_cost_constants.at(2), edit_cost_constants.at(3), edit_cost_constants.at(4), edit_cost_constants.at(5));
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new Constant<UserNodeLabel, UserEdgeLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::CONSTANT. Expected: 6 or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	default:
		throw Error("Selected edit costs unavailable for template parameters UserNodeLabel = " + std::string(typeid(UserNodeLabel).name()) + " and UserEdgeLabel = " + std::string(typeid(UserEdgeLabel).name()) + ".");
	}
	delete_edit_costs_ = true;
}

template<>
void
GEDData<GXLLabel, GXLLabel>::
set_edit_costs_(Options::EditCosts edit_costs, const std::vector<double> & edit_cost_constants) {
	if (delete_edit_costs_) {
		delete edit_costs_;
	}
	switch (edit_costs) {
	case Options::EditCosts::CHEM_1:
		if (edit_cost_constants.size() == 4) {
			edit_costs_ = new CHEM1<GXLLabel, GXLLabel>(edit_cost_constants.at(0), edit_cost_constants.at(1), edit_cost_constants.at(2), edit_cost_constants.at(3)); // @suppress("Symbol is not resolved")
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new CHEM1<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::CHEM_1. Expected: 4 or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	case Options::EditCosts::CHEM_2:
		if (edit_cost_constants.size() == 3) {
			edit_costs_ = new CHEM2<GXLLabel, GXLLabel>(edit_cost_constants.at(0), edit_cost_constants.at(1), edit_cost_constants.at(2)); // @suppress("Symbol is not resolved")
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new CHEM2<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::CHEM_2. Expected: 3 or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	case Options::EditCosts::GREC_1:
		if (edit_cost_constants.size() == 0) {
			edit_costs_ = new GREC1<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::GREC_1. Expected: 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	case Options::EditCosts::GREC_2:
		if (edit_cost_constants.size() == 3) {
			edit_costs_ = new GREC2<GXLLabel, GXLLabel>(edit_cost_constants.at(0), edit_cost_constants.at(1), edit_cost_constants.at(2)); // @suppress("Symbol is not resolved")
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new GREC2<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::GREC_2. Expected: 3 or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	case Options::EditCosts::PROTEIN:
		if (edit_cost_constants.size() == 3) {
			edit_costs_ = new Protein<GXLLabel, GXLLabel>(edit_cost_constants.at(0), edit_cost_constants.at(1), edit_cost_constants.at(2)); // @suppress("Symbol is not resolved")
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new Protein<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::PROTEIN. Expected: 3 or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	case Options::EditCosts::FINGERPRINT:
		if (edit_cost_constants.size() == 3) {
			edit_costs_ = new Fingerprint<GXLLabel, GXLLabel>(edit_cost_constants.at(0), edit_cost_constants.at(1), edit_cost_constants.at(2)); // @suppress("Symbol is not resolved")
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new Fingerprint<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::FINGERPRINT. Expected: 3 or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	case Options::EditCosts::LETTER:
		if (edit_cost_constants.size() == 4) {
			edit_costs_ = new Letter<GXLLabel, GXLLabel>(edit_cost_constants.at(0), edit_cost_constants.at(1), edit_cost_constants.at(2), edit_cost_constants.at(3)); // @suppress("Symbol is not resolved")
		}
		else if (edit_cost_constants.size() == 1) {
			edit_costs_ = new Letter<GXLLabel, GXLLabel>(edit_cost_constants.at(0));
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new Letter<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::LETTER. Expected: 4, 1, or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	case Options::EditCosts::CMU:
		if (edit_cost_constants.size() == 2) {
			edit_costs_ = new CMU<GXLLabel, GXLLabel>(edit_cost_constants.at(0), edit_cost_constants.at(2)); // @suppress("Symbol is not resolved")
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new CMU<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::CMU. Expected: 2 or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	case Options::EditCosts::CONSTANT:
		if (edit_cost_constants.size() == 6) {
			edit_costs_ = new Constant<GXLLabel, GXLLabel>(edit_cost_constants.at(0), edit_cost_constants.at(1), edit_cost_constants.at(2), edit_cost_constants.at(3), edit_cost_constants.at(4), edit_cost_constants.at(5)); // @suppress("Symbol is not resolved")
		}
		else if (edit_cost_constants.size() == 0) {
			edit_costs_ = new Constant<GXLLabel, GXLLabel>();
		}
		else {
			throw Error("Wrong number of constants for selected edit costs ged::Options::EditCosts::CONSTANT. Expected: 6 or 0; actual: " + std::to_string(edit_cost_constants.size()) + ".");
		}
		break;
	}
	delete_edit_costs_ = true;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
set_edit_costs_(EditCosts<UserNodeLabel, UserEdgeLabel> * edit_costs) {
	if (delete_edit_costs_) {
		delete edit_costs_;
	}
	edit_costs_ = edit_costs;
	delete_edit_costs_ = false;
}

template<class UserNodeLabel, class UserEdgeLabel>
LabelID
GEDData<UserNodeLabel, UserEdgeLabel>::
node_label_to_id_(const UserNodeLabel & node_label) {
    auto itr = node_label_ids_.find(node_label);
    if (itr != node_label_ids_.end()) {
        return itr->second;
    }
    node_labels_.push_back(node_label);
    LabelID node_label_id{node_labels_.size()};
    node_label_ids_[node_label] = node_label_id;
    return node_label_id;
}

template<class UserNodeLabel, class UserEdgeLabel>
UserNodeLabel
GEDData<UserNodeLabel, UserEdgeLabel>::
id_to_node_label(LabelID label_id) const {
	if ((label_id > node_labels_.size()) or (label_id == 0)) {
		throw Error("Invalid node label ID " + std::to_string(label_id) + ".");
	}
	return node_labels_.at(label_id - 1);
}


template<class UserNodeLabel, class UserEdgeLabel>
LabelID
GEDData<UserNodeLabel, UserEdgeLabel>::
edge_label_to_id_(const UserEdgeLabel & edge_label) {
    auto itr = edge_label_ids_.find(edge_label);
    if (itr != edge_label_ids_.end()) {
        return itr->second;
    }
    edge_labels_.push_back(edge_label);
    LabelID edge_label_id{edge_labels_.size()};
    edge_label_ids_[edge_label] = edge_label_id;
    return edge_label_id;
}



template<class UserNodeLabel, class UserEdgeLabel>
UserEdgeLabel
GEDData<UserNodeLabel, UserEdgeLabel>::
id_to_edge_label(LabelID label_id) const {
	if ((label_id > edge_labels_.size()) or (label_id == 0)) {
		throw Error("Invalid edge label ID " + std::to_string(label_id) + ".");
	}
	return edge_labels_.at(label_id - 1);
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
GEDData<UserNodeLabel, UserEdgeLabel>::
eager_init_() const {
	return ((init_type_ == Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES) or (init_type_ == Options::InitType::EAGER_WITH_SHUFFLED_COPIES));
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
init_cost_matrices_(bool print_to_stdout) {

	// Update node cost matrix if new node labels have been added to the environment.
	std::size_t size_old_node_costs = node_costs_.num_rows();
	if (size_old_node_costs < node_labels_.size() + 1) {
		DMatrix old_node_costs(node_costs_);
		node_costs_.resize(node_labels_.size() + 1, node_labels_.size() + 1);
		ProgressBar progress(node_labels_.size() + 1);
		if (print_to_stdout) {
			std::cout << "\rInitializing node cost matrix: " << progress << std::flush;
		}
		for (LabelID l_id_lhs = 0; l_id_lhs < node_labels_.size() + 1; l_id_lhs++) {
			for (LabelID l_id_rhs = 0; l_id_rhs < node_labels_.size() + 1; l_id_rhs++) {
				if (l_id_lhs < size_old_node_costs and l_id_rhs < size_old_node_costs) {
					node_costs_(l_id_lhs, l_id_rhs) = old_node_costs(l_id_lhs, l_id_rhs);
				}
				else if (l_id_lhs == l_id_rhs) {
					node_costs_(l_id_lhs, l_id_rhs) = 0.0;
				}
				else if (l_id_lhs == dummy_label()) {
					node_costs_(l_id_lhs, l_id_rhs) = edit_costs_->node_ins_cost_fun(node_labels_.at(l_id_rhs - 1));
				}
				else if (l_id_rhs == dummy_label()) {
					node_costs_(l_id_lhs, l_id_rhs) = edit_costs_->node_del_cost_fun(node_labels_.at(l_id_lhs - 1));
				}
				else {
					node_costs_(l_id_lhs, l_id_rhs) = edit_costs_->node_rel_cost_fun(node_labels_.at(l_id_lhs - 1), node_labels_.at(l_id_rhs - 1));
				}
			}
			if (print_to_stdout) {
				progress.increment();
				std::cout << "\rInitializing node cost matrix: " << progress << std::flush;
			}
		}
		if (print_to_stdout) {
			std::cout << "\n";
		}
	}

	// Update edge cost matrix if new edge labels have been added to the environment.
	std::size_t size_old_edge_costs = edge_costs_.num_rows();
	if (size_old_edge_costs < edge_labels_.size() + 1) {
		DMatrix old_edge_costs(edge_costs_);
		edge_costs_.resize(edge_labels_.size() + 1, edge_labels_.size() + 1);
		ProgressBar progress(edge_labels_.size() + 1);
		if (print_to_stdout) {
			std::cout << "\rInitializing edge cost matrix: " << progress << std::flush;
		}
		for (LabelID l_id_lhs = 0; l_id_lhs < edge_labels_.size() + 1; l_id_lhs++) {
			for (LabelID l_id_rhs = 0; l_id_rhs < edge_labels_.size() + 1; l_id_rhs++) {
				if (l_id_lhs < size_old_edge_costs and l_id_rhs < size_old_edge_costs) {
					edge_costs_(l_id_lhs, l_id_rhs) = old_edge_costs(l_id_lhs, l_id_rhs);
				}
				else if (l_id_lhs == l_id_rhs) {
					edge_costs_(l_id_lhs, l_id_rhs) = 0.0;
				}
				else if (l_id_lhs == dummy_label()) {
					edge_costs_(l_id_lhs, l_id_rhs) = edit_costs_->edge_ins_cost_fun(edge_labels_.at(l_id_rhs - 1));
				}
				else if (l_id_rhs == dummy_label()) {
					edge_costs_(l_id_lhs, l_id_rhs) = edit_costs_->edge_del_cost_fun(edge_labels_.at(l_id_lhs - 1));
				}
				else {
					edge_costs_(l_id_lhs, l_id_rhs) = edit_costs_->edge_rel_cost_fun(edge_labels_.at(l_id_lhs - 1), edge_labels_.at(l_id_rhs - 1));
				}
			}
			if (print_to_stdout) {
				progress.increment();
				std::cout << "\rInitializing edge cost matrix: " << progress << std::flush;
			}
		}
		if (print_to_stdout) {
			std::cout << "\n";
		}
	}

}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDData<UserNodeLabel, UserEdgeLabel>::
num_graphs() const {
	return graphs_.size();
}

template<class UserNodeLabel, class UserEdgeLabel>
const GEDGraph &
GEDData<UserNodeLabel, UserEdgeLabel>::
graph(GEDGraph::GraphID graph_id) const {
	return graphs_.at(graph_id);
}

template<class UserNodeLabel, class UserEdgeLabel>
GEDGraph::NodeID
GEDData<UserNodeLabel, UserEdgeLabel>::
string_to_node_id_(GEDGraph::GraphID graph_id, const std::string & string) const {
	if (string == "DUMMY") {
		return GEDGraph::dummy_node();
	}
	return strings_to_internal_node_ids_.at(graph_id).at(string);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
GEDData<UserNodeLabel, UserEdgeLabel>::
node_id_to_string_(GEDGraph::GraphID graph_id, GEDGraph::NodeID node_id) const {
	if (node_id == GEDGraph::dummy_node()) {
		return "DUMMY";
	}
	return internal_node_ids_to_strings_.at(graph_id).at(node_id);
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
GEDData<UserNodeLabel, UserEdgeLabel>::
shuffled_graph_copies_available() const {
	return ((init_type_ == Options::InitType::EAGER_WITH_SHUFFLED_COPIES) or (init_type_ == Options::InitType::LAZY_WITH_SHUFFLED_COPIES));
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDData<UserNodeLabel, UserEdgeLabel>::
num_graphs_without_shuffled_copies() const {
	return num_graphs_without_shuffled_copies_;
}

template<class UserNodeLabel, class UserEdgeLabel>
GEDGraph::GraphID
GEDData<UserNodeLabel, UserEdgeLabel>::
id_shuffled_graph_copy(GEDGraph::GraphID graph_id) const {
	if (not shuffled_graph_copies_available()) {
		throw Error("No shuffled copy available.");
	}
	if (graph_id >= num_graphs_without_shuffled_copies()) {
		return (graph_id - num_graphs_without_shuffled_copies());
	}
	return (graph_id + num_graphs_without_shuffled_copies());
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
GEDData<UserNodeLabel, UserEdgeLabel>::
is_shuffled_graph_copy(GEDGraph::GraphID graph_id) const {
	return (graph_id >= num_graphs_without_shuffled_copies());
}

template<class UserNodeLabel, class UserEdgeLabel>
std::vector<GEDGraph>::const_iterator
GEDData<UserNodeLabel, UserEdgeLabel>::
begin() const {
	return graphs_.cbegin();
}

template<class UserNodeLabel, class UserEdgeLabel>
std::vector<GEDGraph>::const_iterator
GEDData<UserNodeLabel, UserEdgeLabel>::
end() const {
	return graphs_.cend();
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDData<UserNodeLabel, UserEdgeLabel>::
num_node_labels() const {
	return node_costs_.num_rows();
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDData<UserNodeLabel, UserEdgeLabel>::
num_edge_labels() const {
	return edge_costs_.num_rows();
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDData<UserNodeLabel, UserEdgeLabel>::
max_num_nodes() const {
	return max_num_nodes_;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
GEDData<UserNodeLabel, UserEdgeLabel>::
max_num_edges() const {
	return max_num_edges_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
edge_cost(LabelID label1, LabelID label2) const {
	if (eager_init_()) {
		return edge_costs_(label1, label2);
	}
	if (label1 == label2) {
		return 0.0;
	}
	if (label1 == dummy_label()) {
		return edit_costs_->edge_ins_cost_fun(edge_labels_.at(label2 - 1));
	}
	if (label2 == dummy_label()) {
		return edit_costs_->edge_del_cost_fun(edge_labels_.at(label1 - 1));
	}
	return edit_costs_->edge_rel_cost_fun(edge_labels_.at(label1 - 1), edge_labels_.at(label2 - 1));
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
vectorize_edge_label(LabelID edge_label, std::vector<double> & vector_representation) const {
	edit_costs_->vectorize_edge_label(edge_labels_.at(edge_label - 1), vector_representation);
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
node_cost(LabelID label1, LabelID label2) const {
	if (eager_init_()) {
		return node_costs_(label1, label2);
	}
	if (label1 == label2) {
		return 0.0;
	}
	if (label1 == dummy_label()) {
		return edit_costs_->node_ins_cost_fun(node_labels_.at(label2 - 1));
	}
	if (label2 == dummy_label()) {
		return edit_costs_->node_del_cost_fun(node_labels_.at(label1 - 1));
	}
	return edit_costs_->node_rel_cost_fun(node_labels_.at(label1 - 1), node_labels_.at(label2 - 1));
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
vectorize_node_label(LabelID node_label, std::vector<double> & vector_representation) const {
	edit_costs_->vectorize_node_label(node_labels_.at(node_label - 1), vector_representation);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
save_node_map(const std::string & filename, GEDGraph::NodeID g_id, GEDGraph::NodeID h_id, const NodeMap & node_map, bool append) const {
	std::ofstream ofs;
	if (append) {
		ofs.open(filename, std::ofstream::out | std::ofstream::app);
	}
	else {
		ofs.open(filename, std::ofstream::out);
	}
	ofs << graph_names_.at(g_id) << " " << graph_names_.at(h_id);
	std::vector<NodeMap::Assignment> relation;
	node_map.as_relation(relation);
	for (const auto & assignment : relation) {
		ofs << " " << node_id_to_string_(g_id, assignment.first) << "," << node_id_to_string_(h_id, assignment.second);
	}
	ofs << "\n";
	ofs.close();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
load_node_map(const std::string & filename, GEDGraph::NodeID g_id, GEDGraph::NodeID h_id, NodeMap & node_map) const {
	std::string prefix(graph_names_.at(g_id) + " " + graph_names_.at(h_id) + " ");
	std::ifstream ifs(filename);
	if (not ifs.good()) {
		throw Error("Loading node map from file " + filename + " failed. File cannot be opened.");
	}
	std::string line;
	bool contains_node_map{false};
	while(std::getline(ifs, line)) {
		auto res = std::mismatch(prefix.begin(), prefix.end(), line.begin());
		if (res.first == prefix.end()) {
			contains_node_map = true;
			line = line.substr(prefix.size());
			break;
		}
	}
	ifs.close();
	if (not contains_node_map) {
		throw Error("Loading node map from file " + filename + " failed. File does not contain a node map between the graphs " + graph_names_.at(g_id) + " and " + graph_names_.at(h_id) + ".");
	}
	std::istringstream line_stream(line);
	std::string assignment;
	std::string node_in_g;
	std::string node_in_h;
	while (std::getline(line_stream, assignment, ' ')) {
		std::istringstream assignment_stream(assignment);
		if (not std::getline(assignment_stream, node_in_g, ',')) {
			throw Error("Loading node map from file " + filename + " failed. File has the wrong format. Expected format of lines: \"<graph name 1> <graph name 2> [<original node ID 1>,<original node ID 2>] [...]\"");
		}
		if (not std::getline(assignment_stream, node_in_h, ',')) {
			throw Error("Loading node map from file " + filename + " failed. File has the wrong format. Expected format of lines: \"<graph name 1> <graph name 2> [<original node ID 1>,<original node ID 2>] [...]\"");
		}
		node_map.add_assignment(string_to_node_id_(g_id, node_in_g), string_to_node_id_(h_id, node_in_h));
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
compute_induced_cost(const GEDGraph & g, const GEDGraph & h, NodeMap & node_map) const {
	double cost{0.0};

	// collect node costs
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		cost += node_cost(g.get_node_label(*node), h.get_node_label(node_map.image(*node)));
	}
	for (auto node = h.nodes().first; node != h.nodes().second; node++) {
		GEDGraph::NodeID pre_image{node_map.pre_image(*node)};
		if (pre_image == GEDGraph::dummy_node()) {
			cost += node_cost(g.get_node_label(pre_image), h.get_node_label(*node));
		}
	}

	// collect edge costs
	for (auto ij = g.edges().first; ij != g.edges().second; ij++) {
		cost += edge_cost(g.get_edge_label(*ij), h.get_edge_label(node_map.image(g.tail(*ij)), node_map.image(g.head(*ij))));
	}
	for (auto kl = h.edges().first; kl != h.edges().second; kl++) {
		if (not g.is_edge(node_map.pre_image(h.tail(*kl)), node_map.pre_image(h.head(*kl)))) {
			cost += edge_cost(dummy_label(), h.get_edge_label(*kl));
		}
	}

	node_map.set_induced_cost(cost);
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
swap_cost(const GEDGraph & g, const GEDGraph & h, const NodeMap::Assignment & assignment_1, const NodeMap::Assignment & assignment_2, NodeMap & node_map) const {

	// Get the nodes involved in the swap: {(i,k), (j,l)} -> {(i,l), (j,k)}.
	GEDGraph::NodeID i{assignment_1.first};
	GEDGraph::NodeID k{assignment_1.second};
	GEDGraph::NodeID j{assignment_2.first};
	GEDGraph::NodeID l{assignment_2.second};

	// Collect edges that are incident with i or j in g.
	std::vector<GEDGraph::EdgeID> incident_edges_i_j;
	if (i != GEDGraph::dummy_node()) {
		for (auto edge = g.incident_edges(i).first; edge != g.incident_edges(i).second; edge++) {
			if (g.head(*edge) != j) {
				incident_edges_i_j.push_back(*edge);
			}
		}
	}
	if (j != GEDGraph::dummy_node()) {
		for (auto edge = g.incident_edges(j).first; edge != g.incident_edges(j).second; edge++) {
			if (g.head(*edge) != i) {
				incident_edges_i_j.push_back(*edge);
			}
		}
	}

	// Collect edges that are incident with k or l in h.
	std::vector<GEDGraph::EdgeID> incident_edges_k_l;
	if (k != GEDGraph::dummy_node()) {
		for (auto edge = h.incident_edges(k).first; edge != h.incident_edges(k).second; edge++) {
			if (h.head(*edge) != l) {
				incident_edges_k_l.push_back(*edge);
			}
		}
	}
	if (l != GEDGraph::dummy_node()) {
		for (auto edge = h.incident_edges(l).first; edge != h.incident_edges(l).second; edge++) {
			if (h.head(*edge) != k) {
				incident_edges_k_l.push_back(*edge);
			}
		}
	}

	// Compute swap cost.
	double delta{0.0};

	// Compute node cost delta.
	delta -= node_cost(g.get_node_label(i), h.get_node_label(k));
	delta -= node_cost(g.get_node_label(j), h.get_node_label(l));
	delta += node_cost(g.get_node_label(i), h.get_node_label(l));
	delta += node_cost(g.get_node_label(j), h.get_node_label(k));

	// Compute negative part of edge cost delta.
	for (const auto & edge : incident_edges_i_j) {
		delta -= edge_cost(g.get_edge_label(edge), h.get_edge_label(node_map.image(g.tail(edge)), node_map.image(g.head(edge))));
	}
	for (const auto & edge : incident_edges_k_l) {
		if (not g.is_edge(node_map.pre_image(h.tail(edge)), node_map.pre_image(h.head(edge)))) {
			delta -= edge_cost(dummy_label(), h.get_edge_label(edge));
		}
	}

	// Carry out the swap.
	node_map.add_assignment(i, l);
	node_map.add_assignment(j, k);

	// Compute positive part of edge cost delta.
	for (const auto & edge : incident_edges_i_j) {
		delta += edge_cost(g.get_edge_label(edge), h.get_edge_label(node_map.image(g.tail(edge)), node_map.image(g.head(edge))));
	}
	for (const auto & edge : incident_edges_k_l) {
		if (not g.is_edge(node_map.pre_image(h.tail(edge)), node_map.pre_image(h.head(edge)))) {
			delta += edge_cost(dummy_label(), h.get_edge_label(edge));
		}
	}

	// Undo the swap.
	node_map.add_assignment(i, k);
	node_map.add_assignment(j, l);

	// Return the overall swap cost.
	return delta;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDData<UserNodeLabel, UserEdgeLabel>::
swap_assignments(const NodeMap::Assignment & assignment_1, const NodeMap::Assignment & assignment_2, double swap_cost, NodeMap & node_map) const {

	// Carry out the swap.
	node_map.add_assignment(assignment_1.first, assignment_2.second);
	node_map.add_assignment(assignment_2.first, assignment_1.second);

	// Update the induced cost of the node map.
	node_map.set_induced_cost(node_map.induced_cost() + swap_cost);
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
GEDData<UserNodeLabel, UserEdgeLabel>::
quasimetric_costs() const {
	for (std::size_t label_1{1}; label_1 < num_node_labels(); label_1++) {
		for (std::size_t label_2{1}; label_2 < num_node_labels(); label_2++) {
			if (node_cost(label_1, label_2) > node_cost(label_1, ged::dummy_label()) + node_cost(ged::dummy_label(), label_2)) {
				return false;
			}
		}
	}
	for (std::size_t label_1{1}; label_1 < num_edge_labels(); label_1++) {
		for (std::size_t label_2{1}; label_2 < num_edge_labels(); label_2++) {
			if (edge_cost(label_1, label_2) > edge_cost(label_1, ged::dummy_label()) + edge_cost(ged::dummy_label(), label_2)) {
				return false;
			}
		}
	}
	return true;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
GEDData<UserNodeLabel, UserEdgeLabel>::
quasimetric_costs(const GEDGraph & g, const GEDGraph & h) const {
	std::vector<LabelID> node_labels_g;
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		node_labels_g.push_back(g.get_node_label(*node));
	}
	std::vector<LabelID> node_labels_h;
	for (auto node = h.nodes().first; node != h.nodes().second; node++) {
		node_labels_h.push_back(h.get_node_label(*node));
	}
	std::vector<LabelID> edge_labels_g;
	for (auto edge = g.edges().first; edge != g.edges().second; edge++) {
		edge_labels_g.push_back(g.get_edge_label(*edge));
	}
	std::vector<LabelID> edge_labels_h;
	for (auto edge = h.edges().first; edge != h.edges().second; edge++) {
		edge_labels_h.push_back(h.get_edge_label(*edge));
	}
	for (auto label_1 : node_labels_g) {
		for (auto label_2 : node_labels_h) {
			if (node_cost(label_1, label_2) > node_cost(label_1, ged::dummy_label()) + node_cost(ged::dummy_label(), label_2)) {
				return false;
			}
		}
	}
	for (auto label_1 : edge_labels_g) {
		for (auto label_2 : edge_labels_h) {
			if (edge_cost(label_1, label_2) > edge_cost(label_1, ged::dummy_label()) + edge_cost(ged::dummy_label(), label_2)) {
				return false;
			}
		}
	}
	return true;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_node_edit_cost() const {
	if (eager_init_()) {
		return node_costs_.max();
	}
	double max_cost{std::numeric_limits<double>::min()};
	for (LabelID label_1{0}; label_1 < num_node_labels(); label_1++) {
		for (LabelID label_2{0}; label_2 < num_node_labels(); label_2++) {
			max_cost = std::max(max_cost, node_cost(label_1, label_2));
		}
	}
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_node_del_cost(const GEDGraph & graph) const {
	double max_cost{std::numeric_limits<double>::min()};
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		max_cost = std::max(max_cost, node_cost(graph.get_node_label(*node), dummy_label()));
	}
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_node_del_cost(const GEDGraph & graph) const {
	double min_cost{std::numeric_limits<double>::max()};
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		min_cost = std::min(min_cost, node_cost(graph.get_node_label(*node), dummy_label()));
	}
	return min_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
mean_node_del_cost(const GEDGraph & graph) const {
	if (graph.num_nodes() == 0) {
		return 0.0;
	}
	double mean_cost{0.0};
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		mean_cost += node_cost(graph.get_node_label(*node), dummy_label());
	}
	return mean_cost / static_cast<double>(graph.num_nodes());
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_node_ins_cost(const GEDGraph & graph) const {
	double max_cost{std::numeric_limits<double>::min()};
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		max_cost = std::max(max_cost, node_cost(dummy_label(), graph.get_node_label(*node)));
	}
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_node_ins_cost(const GEDGraph & graph) const {
	double min_cost{std::numeric_limits<double>::max()};
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		min_cost = std::min(min_cost, node_cost(dummy_label(), graph.get_node_label(*node)));
	}
	return min_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
mean_node_ins_cost(const GEDGraph & graph) const {
	if (graph.num_nodes() == 0) {
		return 0.0;
	}
	double mean_cost{0.0};
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		mean_cost += node_cost(dummy_label(), graph.get_node_label(*node));
	}
	return mean_cost / static_cast<double>(graph.num_nodes());
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_node_subs_cost(const GEDGraph & g, const GEDGraph & h) const {
	double max_cost{std::numeric_limits<double>::min()};
	for (auto node_g = g.nodes().first; node_g != g.nodes().second; node_g++) {
		for (auto node_h = h.nodes().first; node_h != h.nodes().second; node_h++) {
			max_cost = std::max(max_cost, node_cost(g.get_node_label(*node_g), h.get_node_label(*node_h)));
		}
	}
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_node_subs_cost(const GEDGraph & g, const GEDGraph & h) const {
	double min_cost{std::numeric_limits<double>::max()};
	for (auto node_g = g.nodes().first; node_g != g.nodes().second; node_g++) {
		for (auto node_h = h.nodes().first; node_h != h.nodes().second; node_h++) {
			if (g.get_node_label(*node_g) != h.get_node_label(*node_h)) {
				min_cost = std::min(min_cost, node_cost(g.get_node_label(*node_g), h.get_node_label(*node_h)));
			}
		}
	}
	return min_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
mean_node_subs_cost(const GEDGraph & g, const GEDGraph & h) const {
	if (g.num_nodes() * h.num_nodes() == 0) {
		return 0.0;
	}
	double mean_cost{0.0};
	for (auto node_g = g.nodes().first; node_g != g.nodes().second; node_g++) {
		for (auto node_h = h.nodes().first; node_h != h.nodes().second; node_h++) {
			mean_cost += node_cost(g.get_node_label(*node_g), h.get_node_label(*node_h));
		}
	}
	return mean_cost / static_cast<double>(g.num_nodes() * h.num_nodes());
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_edge_subs_cost(const GEDGraph & g, const GEDGraph & h) const {
	double max_cost{std::numeric_limits<double>::min()};
	for (auto egde_g = g.edges().first; egde_g != g.edges().second; egde_g++) {
		for (auto edge_h = h.edges().first; edge_h != h.edges().second; edge_h++) {
			max_cost = std::max(max_cost, edge_cost(g.get_edge_label(*egde_g), h.get_edge_label(*edge_h)));
		}
	}
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_edge_subs_cost(const GEDGraph & g, const GEDGraph & h) const {
	double min_cost{std::numeric_limits<double>::max()};
	for (auto egde_g = g.edges().first; egde_g != g.edges().second; egde_g++) {
		for (auto edge_h = h.edges().first; edge_h != h.edges().second; edge_h++) {
			if (g.get_edge_label(*egde_g) != h.get_edge_label(*edge_h)) {
				min_cost = std::min(min_cost, edge_cost(g.get_edge_label(*egde_g), h.get_edge_label(*edge_h)));
			}
		}
	}
	return min_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
mean_edge_subs_cost(const GEDGraph & g, const GEDGraph & h) const {
	if (g.num_edges() * h.num_edges() == 0) {
		return 0.0;
	}
	double mean_cost{0.0};
	for (auto edge_g = g.edges().first; edge_g != g.edges().second; edge_g++) {
		for (auto edge_h = h.edges().first; edge_h != h.edges().second; edge_h++) {
			mean_cost += edge_cost(g.get_edge_label(*edge_g), h.get_edge_label(*edge_h));
		}
	}
	return mean_cost / static_cast<double>(g.num_edges() * h.num_edges());
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_edge_edit_cost() const {
	if (eager_init_()) {
		return edge_costs_.max();
	}
	double max_cost{std::numeric_limits<double>::min()};
	for (LabelID label_1{0}; label_1 < num_edge_labels(); label_1++) {
		for (LabelID label_2{0}; label_2 < num_edge_labels(); label_2++) {
			max_cost = std::max(max_cost, edge_cost(label_1, label_2));
		}
	}
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_edge_del_cost(const GEDGraph & graph) const {
	double max_cost{std::numeric_limits<double>::min()};
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		max_cost = std::max(max_cost, edge_cost(graph.get_edge_label(*edge), dummy_label()));
	}
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
mean_edge_del_cost(const GEDGraph & graph) const {
	if (graph.num_edges() == 0) {
		return 0.0;
	}
	double mean_cost{0.0};
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		mean_cost += edge_cost(graph.get_edge_label(*edge), dummy_label());
	}
	return mean_cost / static_cast<double>(graph.num_edges());
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_edge_del_cost(const GEDGraph & graph) const {
	double min_cost{std::numeric_limits<double>::max()};
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		min_cost = std::min(min_cost, edge_cost(graph.get_edge_label(*edge), dummy_label()));
	}
	return min_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_edge_ins_cost(const GEDGraph & graph) const {
	double max_cost{std::numeric_limits<double>::min()};
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		max_cost = std::max(max_cost, edge_cost(dummy_label(), graph.get_edge_label(*edge)));
	}
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_edge_ins_cost(const GEDGraph & graph) const {
	double min_cost{std::numeric_limits<double>::max()};
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		min_cost = std::min(min_cost, edge_cost(dummy_label(), graph.get_edge_label(*edge)));
	}
	return min_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
mean_edge_ins_cost(const GEDGraph & graph) const {
	if (graph.num_edges() == 0) {
		return 0.0;
	}
	double mean_cost{0.0};
	for (auto edge = graph.edges().first; edge != graph.edges().second; edge++) {
		mean_cost += edge_cost(dummy_label(), graph.get_edge_label(*edge));
	}
	return mean_cost / static_cast<double>(graph.num_edges());
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_edit_cost(const GEDGraph & g, const GEDGraph & h) const {
	double max_cost{max_edge_del_cost(g)};
	max_cost = std::max(max_cost, max_edge_ins_cost(h));
	max_cost = std::max(max_cost, max_edge_subs_cost(g, h));
	max_cost = std::max(max_cost, max_node_del_cost(g));
	max_cost = std::max(max_cost, max_node_ins_cost(h));
	max_cost = std::max(max_cost, max_node_subs_cost(g, h));
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_edit_cost(const GEDGraph & g, const GEDGraph & h) const {
	double min_cost{min_edge_del_cost(g)};
	min_cost = std::min(min_cost, min_edge_ins_cost(h));
	min_cost = std::min(min_cost, min_edge_subs_cost(g, h));
	min_cost = std::min(min_cost, min_node_del_cost(g));
	min_cost = std::min(min_cost, min_node_ins_cost(h));
	min_cost = std::min(min_cost, min_node_subs_cost(g, h));
	return min_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_node_edit_cost(const GEDGraph & g, const GEDGraph & h) const {
	double max_cost{max_node_del_cost(g)};
	max_cost = std::max(max_cost, max_node_ins_cost(h));
	max_cost = std::max(max_cost, max_node_subs_cost(g, h));
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_node_edit_cost(const GEDGraph & g, const GEDGraph & h) const {
	double min_cost{min_node_del_cost(g)};
	min_cost = std::min(min_cost, min_node_ins_cost(h));
	min_cost = std::min(min_cost, min_node_subs_cost(g, h));
	return min_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
max_edge_edit_cost(const GEDGraph & g, const GEDGraph & h) const {
	double max_cost{max_edge_del_cost(g)};
	max_cost = std::max(max_cost, max_edge_ins_cost(h));
	max_cost = std::max(max_cost, max_edge_subs_cost(g, h));
	return max_cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDData<UserNodeLabel, UserEdgeLabel>::
min_edge_edit_cost(const GEDGraph & g, const GEDGraph & h) const {
	double min_cost{min_edge_del_cost(g)};
	min_cost = std::min(min_cost, min_edge_ins_cost(h));
	min_cost = std::min(min_cost, min_edge_subs_cost(g, h));
	return min_cost;
}

}

#endif /* SRC_ENV_GED_DATA_IPP_ */
