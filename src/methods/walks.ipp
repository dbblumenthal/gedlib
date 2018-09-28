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
 * @file  walks.ipp
 * @brief ged::Walks class definition
 */

#ifndef SRC_METHODS_WALKS_IPP_
#define SRC_METHODS_WALKS_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
Walks<UserNodeLabel, UserEdgeLabel>::
~Walks() {}

template<class UserNodeLabel, class UserEdgeLabel>
Walks<UserNodeLabel, UserEdgeLabel>::
Walks(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
depth_{3},
min_depth_{1},
max_depth_{5},
infile_(""),
outfile_(""),
adj_graphs_() {
	this->compute_lower_bound_ = false;
}

// Definitions of member functions inherited from LSAPEBasedMethod.
template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {
	depth_ = 3;
	min_depth_ = 1;
	max_depth_ = 5;
	infile_ = std::string("");
	outfile_ = std::string("");
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Walks<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	return "[--depth-range <arg>] [--load <arg>] [--save <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Walks<UserNodeLabel, UserEdgeLabel>::
lsape_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "load") {
		infile_ = arg;
		return true;
	}
	else if (option == "save") {
		outfile_ = arg;
		return true;
	}
	else if (option == "depth-range") {
		std::stringstream depth_range(arg);
		std::string min_depth, max_depth;
		if (std::getline(depth_range, min_depth, ',') and std::getline(depth_range, max_depth, ',')) {
			try {
				min_depth_ = std::stoi(min_depth);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option depth-range. Usage: options = \"[--depth-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
			try {
				max_depth_ = std::stoi(max_depth);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option depth-range. Usage: options = \"[--depth-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
			if ((min_depth_ > max_depth_) or (min_depth_ <= 0) or (max_depth_ <= 0)) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option depth-range. Usage: options = \"[--depth-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
			}
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option depth-range. Usage: options = \"[--depth-range <smaller convertible to int greater 0>,<larger convertible to int greater 0>] [...]");
		}
		return true;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {
	adj_graphs_[graph.id()] = AdjGraph_(graph);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
lsape_pre_graph_init_(bool called_at_runtime) {
	if (load_config_file_()) {
		std::map<std::string, std::string> options;
		util::parse_config_file(infile_, options);
		depth_ = std::stod(options.at("depth"));
	}
	else {
		depth_ = (min_depth_ + max_depth_) / 2;
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
lsape_init_() {

	// Return, if a configuration file is provided or if the minimal depth equals the maximal depth.
	if (load_config_file_() or (min_depth_ == max_depth_)) {
		return;
	}

	// Find the best depth.
	std::size_t num_runs{(max_depth_ - min_depth_ + 1) * (this->ged_data_.num_graphs() * this->ged_data_.num_graphs())};
	ProgressBar progress_bar(num_runs);
	std::cout << "\r" << progress_bar << std::flush;
	double best_avg_ub{std::numeric_limits<double>::infinity()};
	std::size_t best_depth{min_depth_};
	LSAPESolver lsape_solver;
	lsape_solver.set_model(this->lsape_model_);
	for (depth_ = min_depth_; depth_ <= max_depth_; depth_++) {
		double avg_ub{0.0};
		for (auto g = this->ged_data_.begin(); g != this->ged_data_.end(); g++) {
			if (this->ged_data_.is_shuffled_graph_copy(g->id())) {
				continue;
			}
			for (auto h = this->ged_data_.begin(); h != this->ged_data_.end(); h++) {
				if (this->ged_data_.is_shuffled_graph_copy(h->id())) {
					continue;
				}
				NodeMap node_map(g->num_nodes(), h->num_nodes());
				DMatrix lsape_instance(g->num_nodes() + 1, h->num_nodes() + 1, 0.0);
				lsape_solver.set_problem(&lsape_instance);
				if (this->ged_data_.shuffled_graph_copies_available() and (g->id() == h->id())) {
					GEDGraph::GraphID id_shuffled_graph_copy{this->ged_data_.id_shuffled_graph_copy(h->id())};
					lsape_populate_instance_(*g, this->ged_data_.graph(id_shuffled_graph_copy), lsape_instance);
					lsape_solver.solve();
					util::construct_node_map_from_solver(lsape_solver, node_map);
					this->ged_data_.compute_induced_cost(*g, this->ged_data_.graph(id_shuffled_graph_copy), node_map);
				}
				else {
					lsape_populate_instance_(*g, *h, lsape_instance);
					lsape_solver.solve();
					util::construct_node_map_from_solver(lsape_solver, node_map);
					this->ged_data_.compute_induced_cost(*g, *h, node_map);
				}
				avg_ub += node_map.induced_cost();
				progress_bar.increment();
				std::cout << "\r" << progress_bar << std::flush;
			}
		}
		if (avg_ub < best_avg_ub) {
			best_avg_ub = avg_ub;
			best_depth = depth_;
		}
	}
	std::cout << "\n";
	depth_ = best_depth;

	// Save the found depth.
	if (outfile_ != "") {
		std::map<std::string, std::string> options;
		options["depth"] = std::to_string(depth_);
		util::save_as_config_file(outfile_, options);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {

	// Initialize adjacency graphs and product graph.
	std::set<LabelID> node_labels;
	init_node_labels_(g, h, node_labels);
	AdjGraph_ adj_graph_g(adj_graphs_.at(g.id()));
	AdjGraph_ adj_graph_h(adj_graphs_.at(h.id()));
	ProductGraph_ product_graph(g, h);

	// Compute number of unmatched walks.
	adj_graph_g.compute_num_walks_(node_labels, depth_);
	adj_graph_h.compute_num_walks_(node_labels, depth_);
	product_graph.compute_num_unmatched_walks_(node_labels, adj_graph_g, adj_graph_h, depth_);

	// Collect mean costs.
	double mean_node_subs_cost{this->ged_data_.mean_node_subs_cost(g, h)};
	double mean_edge_subs_cost{this->ged_data_.mean_edge_subs_cost(g, h)};
	double mean_node_ins_del_cost{this->ged_data_.mean_node_ins_cost(h) * static_cast<double>(h.num_nodes())};
	mean_node_ins_del_cost += this->ged_data_.mean_node_del_cost(g) * static_cast<double>(g.num_nodes());
	if (g.num_nodes() + h.num_nodes() > 0) {
		mean_node_ins_del_cost /= static_cast<double>(g.num_nodes() + h.num_nodes());
	}
	double mean_edge_ins_del_cost{this->ged_data_.mean_edge_ins_cost(h) * static_cast<double>(h.num_edges())};
	mean_edge_ins_del_cost += this->ged_data_.mean_edge_del_cost(g) * static_cast<double>(g.num_edges());
	if (g.num_edges() + h.num_edges() > 0) {
		mean_edge_ins_del_cost /= static_cast<double>(g.num_edges() + h.num_edges());
	}

	// Populate master problem.
	double depth{static_cast<double>(depth_)};
#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row = 0; row < master_problem.num_rows(); row++) {
		for (std::size_t col = 0; col < master_problem.num_cols(); col++) {
			if ((row < g.num_nodes()) and (col < h.num_nodes())) {
				double delta_1{g.get_node_label(adj_graph_g.node(row)) == h.get_node_label(adj_graph_h.node(col)) ? 0.0 : 1.0};
				master_problem(row,col) = product_graph.num_substituted_walks_ending_at_same_label(row, col) * ((delta_1 + depth - 1) * mean_node_subs_cost + depth * mean_edge_subs_cost);
				master_problem(row,col) += product_graph.num_substituted_walks_ending_at_different_labels(row, col) * ((delta_1 + depth) * mean_node_subs_cost + depth * mean_edge_subs_cost);
				master_problem(row,col) += product_graph.num_inserted_or_deleted_walks(row, col) * ((delta_1 + depth) * mean_node_ins_del_cost + depth * mean_edge_ins_del_cost);
			}
			else if (row < g.num_nodes()) {
				master_problem(row, h.num_nodes()) = product_graph.num_inserted_or_deleted_walks(row, master_problem.num_cols() - 1) * ((1 + depth) * mean_node_ins_del_cost + depth * mean_edge_ins_del_cost);
			}
			else if (col < h.num_nodes()) {
				master_problem(g.num_nodes(), col) = product_graph.num_inserted_or_deleted_walks(master_problem.num_rows() - 1, col) * ((1 + depth) * mean_node_ins_del_cost + depth * mean_edge_ins_del_cost);
			}
		}
	}
}

// === Definition of private class AdjGraph_. ===
template<class UserNodeLabel, class UserEdgeLabel>
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
AdjGraph_(const AdjGraph_ & adj_graph) :
adj_matrix_(adj_graph.adj_matrix_),
num_walks_from_nodes_to_labels_(adj_graph.num_walks_from_nodes_to_labels_),
nodes_(adj_graph.nodes_),
inverse_label_index_(adj_graph.inverse_label_index_){}

template<class UserNodeLabel, class UserEdgeLabel>
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
AdjGraph_() :
adj_matrix_(),
num_walks_from_nodes_to_labels_(),
nodes_(),
inverse_label_index_() {}

template<class UserNodeLabel, class UserEdgeLabel>
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
AdjGraph_(const GEDGraph & g) :
adj_matrix_(g.num_nodes(), g.num_nodes()),
num_walks_from_nodes_to_labels_(),
nodes_(),
inverse_label_index_() {
	std::size_t id{0};
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		nodes_[id] = *node;
		LabelID label{g.get_node_label(*node)};
		if (inverse_label_index_.find(label) == inverse_label_index_.end()) {
			inverse_label_index_[label] = std::vector<std::size_t>{id++};
		}
		else {
			inverse_label_index_[label].push_back(id++);
		}
	}
	for (std::size_t row{0}; row < adj_matrix_.num_rows(); row++) {
		for (std::size_t col{0}; col < adj_matrix_.num_cols(); col++) {
			if (not g.is_edge(node(row), node(col))) {
				adj_matrix_(row, col) = 0;
				continue;
			}
			adj_matrix_(row, col) = 1;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
operator=(const AdjGraph_ & adj_graph) {
	adj_matrix_ = adj_graph.adj_matrix_;
	num_walks_from_nodes_to_labels_ = adj_graph.num_walks_from_nodes_to_labels_;
	nodes_ = adj_graph.nodes_;
	inverse_label_index_ = adj_graph.inverse_label_index_;
}

template<class UserNodeLabel, class UserEdgeLabel>
GEDGraph::NodeID
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
node(std::size_t node_id) const {
	return nodes_.at(node_id);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::pair<std::vector<std::size_t>::const_iterator, std::vector<std::size_t>::const_iterator>
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
nodes_with_label(LabelID label) const {
	if (inverse_label_index_.find(label) != inverse_label_index_.end()) {
		return std::make_pair(inverse_label_index_.at(label).cbegin(), inverse_label_index_.at(label).cend());
	}
	return std::make_pair(inverse_label_index_.begin()->second.cend(), inverse_label_index_.begin()->second.cend());
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
size() const {
	return nodes_.size();
}

template<class UserNodeLabel, class UserEdgeLabel>
const typename Walks<UserNodeLabel, UserEdgeLabel>::Histogram_ &
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
num_walks_from_node_to_labels(std::size_t node_id) const {
	return num_walks_from_nodes_to_labels_.at(node_id);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
compute_num_walks_(const std::set<LabelID> & node_labels, std::size_t depth) {
	adj_matrix_.power(depth);
	for (std::size_t start_node{0}; start_node < size(); start_node++) {
		Histogram_ num_walks_from_node_to_labels;
		for (auto label = node_labels.begin(); label != node_labels.end(); label++) {
			num_walks_from_node_to_labels[*label] = 0.0;
			for (auto end_node = nodes_with_label(*label).first; end_node != nodes_with_label(*label).second; end_node++) {
				num_walks_from_node_to_labels[*label] += adj_matrix_(start_node, *end_node);
			}
		}
		num_walks_from_nodes_to_labels_.push_back(num_walks_from_node_to_labels);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Walks<UserNodeLabel, UserEdgeLabel>::
AdjGraph_ ::
operator() (std::size_t row, std::size_t col) const {
	return static_cast<double>(adj_matrix_(row, col));
}

// === Definition of private class ProductGraph_. ===
template<class UserNodeLabel, class UserEdgeLabel>
Walks<UserNodeLabel, UserEdgeLabel>::
ProductGraph_ ::
ProductGraph_(const GEDGraph & g, const GEDGraph & h) :
adj_matrix_(),
num_substituted_walks_ending_at_same_label_(g.num_nodes(), h.num_nodes(), 0.0),
num_substituted_walks_ending_at_different_labels_(g.num_nodes(), h.num_nodes(), 0.0),
num_inserted_or_deleted_walks_(g.num_nodes() + 1, h.num_nodes() + 1, 0.0),
nodes_(),
inverse_label_index_(),
product_node_to_id_() {
	std::size_t id{0};
	for (auto node_g = g.nodes().first; node_g != g.nodes().second; node_g++) {
		LabelID label_g{g.get_node_label(*node_g)};
		for (auto node_h = h.nodes().first; node_h != h.nodes().second; node_h++) {
			LabelID label_h{h.get_node_label(*node_h)};
			if (label_g == label_h) {
				nodes_.push_back(std::make_pair(*node_g, *node_h));
				if (inverse_label_index_.find(label_g) == inverse_label_index_.end()) {
					inverse_label_index_[label_g] = std::vector<std::size_t>{id++};
				}
				else {
					inverse_label_index_[label_g].push_back(id++);
				}
			}
		}
	}
	adj_matrix_.resize(nodes_.size(), nodes_.size());
	for (std::size_t row{0}; row < adj_matrix_.num_rows(); row++) {
		for (std::size_t col{0}; col < adj_matrix_.num_cols(); col++) {
			if (not g.is_edge(nodes_.at(row).first, nodes_.at(col).first)) {
				adj_matrix_(row, col) = 0;
				continue;
			}
			if (not h.is_edge(nodes_.at(row).second, nodes_.at(col).second)) {
				adj_matrix_(row, col) = 0;
				continue;
			}
			if (g.get_edge_label(g.get_edge(nodes_.at(row).first, nodes_.at(col).first)) != h.get_edge_label(h.get_edge(nodes_.at(row).second, nodes_.at(col).second))) {
				adj_matrix_(row, col) = 0;
				continue;
			}
			adj_matrix_(row, col) = 1;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
std::pair<std::vector<std::size_t>::const_iterator, std::vector<std::size_t>::const_iterator>
Walks<UserNodeLabel, UserEdgeLabel>::
ProductGraph_ ::
nodes_with_label(LabelID label) const {
	if (inverse_label_index_.find(label) != inverse_label_index_.end()) {
		return std::make_pair(inverse_label_index_.at(label).cbegin(), inverse_label_index_.at(label).cend());
	}
	return std::make_pair(inverse_label_index_.begin()->second.cend(), inverse_label_index_.begin()->second.cend());
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
Walks<UserNodeLabel, UserEdgeLabel>::
ProductGraph_ ::
node_id(GEDGraph::NodeID node_g, GEDGraph::NodeID node_h) const {
	ProductNode_ node{std::make_pair(node_g, node_h)};
	if (product_node_to_id_.find(node) != product_node_to_id_.end()) {
		return product_node_to_id_.at(node);
	}
	return undefined_();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
ProductGraph_ ::
compute_num_unmatched_walks_(const std::set<LabelID> & node_labels, const AdjGraph_ & adj_graph_g, const AdjGraph_ & adj_graph_h, std::size_t depth) {
	adj_matrix_.power(depth);
	for (std::size_t start_node_g{0}; start_node_g < adj_graph_g.size(); start_node_g++) {
		for (std::size_t start_node_h{0}; start_node_h < adj_graph_h.size(); start_node_h++) {

			// compute histograms for product graph.
			Histogram_ num_unmatched_walks_from_node_to_labels_g(adj_graph_g.num_walks_from_node_to_labels(start_node_g));
			Histogram_ num_unmatched_walks_from_node_to_labels_h(adj_graph_h.num_walks_from_node_to_labels(start_node_h));
			std::size_t start_node_product{node_id(adj_graph_g.node(start_node_g), adj_graph_h.node(start_node_h))};
			if (start_node_product != undefined_()) {
				for (auto label = node_labels.begin(); label != node_labels.end(); label++) {
					double num_shared_paths_from_node_to_label{0.0};
					for (auto end_node_product = nodes_with_label(*label).first; end_node_product != nodes_with_label(*label).second; end_node_product++) {
						num_shared_paths_from_node_to_label += static_cast<double>(adj_matrix_(start_node_product, *end_node_product));
					}
					num_shared_paths_from_node_to_label = std::sqrt(num_shared_paths_from_node_to_label);
					num_shared_paths_from_node_to_label = std::min(num_shared_paths_from_node_to_label, num_unmatched_walks_from_node_to_labels_g.at(*label));
					num_shared_paths_from_node_to_label = std::min(num_shared_paths_from_node_to_label, num_unmatched_walks_from_node_to_labels_h.at(*label));
					num_unmatched_walks_from_node_to_labels_g[*label] -= num_shared_paths_from_node_to_label;
					num_unmatched_walks_from_node_to_labels_h[*label] -= num_shared_paths_from_node_to_label;
				}
			}

			// compute subsitution costs matrices.
			double subs_same_sum{0.0};
			double r_row_col{0.0};
			double r_col_row{0.0};
			for (auto label = node_labels.begin(); label != node_labels.end(); label++) {
				double subs_same{std::min(num_unmatched_walks_from_node_to_labels_g.at(*label), num_unmatched_walks_from_node_to_labels_h.at(*label))};
				subs_same_sum += subs_same;
				r_row_col += num_unmatched_walks_from_node_to_labels_g.at(*label) - subs_same;
				r_col_row += num_unmatched_walks_from_node_to_labels_h.at(*label) - subs_same;
			}
			num_substituted_walks_ending_at_same_label_(start_node_g, start_node_h) = subs_same_sum;
			num_substituted_walks_ending_at_different_labels_(start_node_g, start_node_h) = std::min(r_row_col, r_col_row);
			num_inserted_or_deleted_walks_(start_node_g, start_node_h) = std::fabs(r_row_col - r_col_row);
		}
	}

	// compute deletion costs matrix.
	for (std::size_t node_id_1{0}; node_id_1 < adj_graph_g.size(); node_id_1++) {
		for (std::size_t node_id_2{0}; node_id_2 < adj_graph_g.size(); node_id_2++) {
			num_inserted_or_deleted_walks_(node_id_1, adj_graph_h.size()) += adj_graph_g(node_id_1, node_id_2);
		}
	}

	// compute insertion cost matrix.
	for (std::size_t node_id_1{0}; node_id_1 < adj_graph_h.size(); node_id_1++) {
		for (std::size_t node_id_2{0}; node_id_2 < adj_graph_h.size(); node_id_2++) {
			num_inserted_or_deleted_walks_(adj_graph_g.size(), node_id_1) += adj_graph_h(node_id_1, node_id_2);
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Walks<UserNodeLabel, UserEdgeLabel>::
ProductGraph_ ::
num_substituted_walks_ending_at_same_label(std::size_t row, std::size_t col) const {
	return num_substituted_walks_ending_at_same_label_(row, col);
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Walks<UserNodeLabel, UserEdgeLabel>::
ProductGraph_ ::
num_substituted_walks_ending_at_different_labels(std::size_t row, std::size_t col) const {
	return num_substituted_walks_ending_at_different_labels_(row, col);
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Walks<UserNodeLabel, UserEdgeLabel>::
ProductGraph_ ::
num_inserted_or_deleted_walks(std::size_t row, std::size_t col) const {
	return num_inserted_or_deleted_walks_(row, col);
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Walks<UserNodeLabel, UserEdgeLabel>::
init_node_labels_(const GEDGraph & g, const GEDGraph & h, std::set<LabelID> & node_labels) const {
	node_labels.clear();
	for (auto node_g = g.nodes().first; node_g != g.nodes().second; node_g++) {
		node_labels.insert(g.get_node_label(*node_g));
	}
	for (auto node_h = h.nodes().first; node_h != h.nodes().second; node_h++) {
		node_labels.insert(h.get_node_label(*node_h));
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Walks<UserNodeLabel, UserEdgeLabel>::
load_config_file_() const {
	return (infile_ != "");
}

}

#endif /* SRC_METHODS_WALKS_IPP_ */
