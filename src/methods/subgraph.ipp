/*!
 * @file subgraph.ipp
 * @brief ged::Subgraph class definition.
 */

#ifndef SRC_METHODS_SUBGRAPH_IPP_
#define SRC_METHODS_SUBGRAPH_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===

template<class UserNodeLabel, class UserEdgeLabel>
Subgraph<UserNodeLabel, UserEdgeLabel>::
~Subgraph() {}

template<class UserNodeLabel, class UserEdgeLabel>
Subgraph<UserNodeLabel, UserEdgeLabel>::
Subgraph(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
depth_{2},
min_depth_{1},
max_depth_{5},
infile_(""),
outfile_(""),
exact_options_("--search-method DFS --map-root-to-root TRUE --time-limit 0.001"),
subgraphs_() {}

// === Definitions of member functions inherited from LSAPEBasedMethod.
template<class UserNodeLabel, class UserEdgeLabel>
void
Subgraph<UserNodeLabel, UserEdgeLabel>::
lsape_set_default_options_() {
	depth_ = 2;
	min_depth_ = 1;
	max_depth_ = 5;
	infile_ = std::string("");
	outfile_ = std::string("");
	exact_options_ = std::string("--search-method DFS --map-root-to-root TRUE --time-limit 0.001");
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Subgraph<UserNodeLabel, UserEdgeLabel>::
lsape_valid_options_string_() const {
	return "[--depth-range <arg>] [--load <arg>] [--save <arg>] [--time-limit-subproblem <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Subgraph<UserNodeLabel, UserEdgeLabel>::
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
	else if (option == "time-limit-subproblem") {
		try {
			std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option time-limit-subproblem. Usage: options = \"[--time-limit-subproblem <convertible to double>] [...]");
		}
		exact_options_ = std::string("--search-method DFS --map-root-to-root TRUE --time-limit ") + arg;
		return true;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Subgraph<UserNodeLabel, UserEdgeLabel>::
lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) {

#ifdef _OPENMP
	omp_set_num_threads(this->num_threads_ - 1);
#pragma omp parallel for if(this->num_threads_ > 1)
#endif
	for (std::size_t row_in_master = 0; row_in_master < master_problem.num_rows(); row_in_master++) {
		Exact<UserNodeLabel, UserEdgeLabel> exact(this->ged_data_);
		exact.set_options(exact_options_);
		for (std::size_t col_in_master = 0; col_in_master < master_problem.num_cols(); col_in_master++) {
			if ((row_in_master < g.num_nodes()) and (col_in_master < h.num_nodes())) {
				master_problem(row_in_master, col_in_master) = compute_substitution_cost_(g, h, row_in_master, col_in_master, exact);
			}
			else if (row_in_master < g.num_nodes()) {
				master_problem(row_in_master, h.num_nodes()) = compute_deletion_cost_(g, row_in_master);
			}
			else if (col_in_master < h.num_nodes()) {
				master_problem(g.num_nodes(), col_in_master) = compute_insertion_cost_(h, col_in_master);
			}
		}
	}

}

template<class UserNodeLabel, class UserEdgeLabel>
void
Subgraph<UserNodeLabel, UserEdgeLabel>::
lsape_init_graph_(const GEDGraph & graph) {
	build_subgraphs_(graph);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Subgraph<UserNodeLabel, UserEdgeLabel>::
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
Subgraph<UserNodeLabel, UserEdgeLabel>::
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
		for (auto graph = this->ged_data_.begin(); graph != this->ged_data_.end(); graph++) {
			build_subgraphs_(*graph);
		}
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
	for (auto graph = this->ged_data_.begin(); graph != this->ged_data_.end(); graph++) {
		build_subgraphs_(*graph);
	}

	// Save the found depth.
	if (outfile_ != "") {
		std::map<std::string, std::string> options;
		options["depth"] = std::to_string(depth_);
		util::save_as_config_file(outfile_, options);
	}
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
bool
Subgraph<UserNodeLabel, UserEdgeLabel>::
load_config_file_() const {
	return (infile_ != "");
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Subgraph<UserNodeLabel, UserEdgeLabel>::
compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k, Exact<UserNodeLabel, UserEdgeLabel> & exact) const {
	const GEDGraph & subgraph_i{subgraphs_.at(subgraph_id_(g, i))};
	const GEDGraph & subgraph_k{subgraphs_.at(subgraph_id_(h, k))};
	Result result;
	exact.run_as_util(subgraph_i, subgraph_k, result);
	return result.upper_bound();
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Subgraph<UserNodeLabel, UserEdgeLabel>::
compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const {
	const GEDGraph & subgraph_i{subgraphs_.at(subgraph_id_(g, i))};
	GEDGraph dummy_graph;
	NodeMap matching(subgraph_i.num_nodes(), 0);
	for (auto node = subgraph_i.nodes().first; node != subgraph_i.nodes().second; node++) {
		matching.add_assignment(*node, GEDGraph::dummy_node());
	}
	this->ged_data_.compute_induced_cost(subgraph_i, dummy_graph, matching);
	return matching.induced_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Subgraph<UserNodeLabel, UserEdgeLabel>::
compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const {
	const GEDGraph & subgraph_k{subgraphs_.at(subgraph_id_(h, k))};
	GEDGraph dummy_graph;
	NodeMap matching(0, subgraph_k.num_nodes());
	for (auto node = subgraph_k.nodes().first; node != subgraph_k.nodes().second; node++) {
		matching.add_assignment(GEDGraph::dummy_node(), *node);
	}
	this->ged_data_.compute_induced_cost(dummy_graph, subgraph_k, matching);
	return matching.induced_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Subgraph<UserNodeLabel, UserEdgeLabel>::
build_subgraphs_(const GEDGraph & graph) {
	for (auto node = graph.nodes().first; node != graph.nodes().second; node++) {
		build_subgraph_(graph, *node);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Subgraph<UserNodeLabel, UserEdgeLabel>::
build_subgraph_(const GEDGraph & graph, GEDGraph::NodeID root_node) {
	GEDGraph::GraphID subgraph_id{subgraph_id_(graph, root_node)};
	subgraphs_[subgraph_id] = GEDGraph(subgraph_id);
	GEDGraph::NodeID root_in_subgraph{subgraphs_.at(subgraph_id).add_node()};
	subgraphs_.at(subgraph_id).set_label(root_in_subgraph, graph.get_node_label(root_node));
	GEDGraph::NodeNodeMap ids_in_subgraph;
	ids_in_subgraph[root_node] = root_in_subgraph;
	build_subgraph_dfs_(graph, root_node, 0, ids_in_subgraph, subgraphs_.at(subgraph_id));
	subgraphs_.at(subgraph_id).setup_adjacency_matrix();
}

template<class UserNodeLabel, class UserEdgeLabel>
GEDGraph::GraphID
Subgraph<UserNodeLabel, UserEdgeLabel>::
subgraph_id_(const GEDGraph & graph, GEDGraph::NodeID node) const {
	return (this->ged_data_.max_num_nodes() * graph.id()) + node;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Subgraph<UserNodeLabel, UserEdgeLabel>::
build_subgraph_dfs_(const GEDGraph & graph, GEDGraph::NodeID current_node, std::size_t depth_current_node, GEDGraph::NodeNodeMap & ids_in_subgraph, GEDGraph & subgraph) {
	if (depth_current_node >= depth_) {
		return;
	}
	GEDGraph::NodeID current_node_in_subgraph{ids_in_subgraph.at(current_node)};
	for (auto edge = graph.incident_edges(current_node).first; edge != graph.incident_edges(current_node).second; edge++) {
		GEDGraph::NodeID next_node{graph.head(*edge)};
		GEDGraph::NodeID next_node_in_subgraph;
		bool found_new_node{false};
		if (ids_in_subgraph.find(next_node) == ids_in_subgraph.end()) {
			found_new_node = true;
			next_node_in_subgraph = subgraph.add_node();
			subgraph.set_label(next_node_in_subgraph, graph.get_node_label(next_node));
			ids_in_subgraph[next_node] = next_node_in_subgraph;
		}
		else {
			next_node_in_subgraph = ids_in_subgraph.at(next_node);
		}
		if (not subgraph.safe_is_edge(current_node_in_subgraph, next_node_in_subgraph)) {
			GEDGraph::EdgeID next_edge{subgraph.add_edge(current_node_in_subgraph, next_node_in_subgraph)};
			subgraph.set_label(next_edge, graph.get_edge_label(*edge));
		}
		if (found_new_node) {
			build_subgraph_dfs_(graph, next_node, depth_current_node + 1, ids_in_subgraph, subgraph);
		}
	}
}

}

#endif /* SRC_METHODS_SUBGRAPH_IPP_ */
