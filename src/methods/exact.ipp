/*!
 * @file  exact.ipp
 * @brief Exact class definition.
 */

#ifndef SRC_METHODS_EXACT_IPP_
#define SRC_METHODS_EXACT_IPP_

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
Exact<UserNodeLabel, UserEdgeLabel>::
~Exact() {}

template<class UserNodeLabel, class UserEdgeLabel>
Exact<UserNodeLabel, UserEdgeLabel>::
Exact(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
lsape_model_{LSAPESolver::ECBP},
search_method_{DFS},
lower_bound_method_{BRANCH_FAST},
num_threads_{1},
time_limit_in_sec_{0.0},
map_root_to_root_{false},
sorted_edges_(),
best_feasible_(this),
open_(),
omega_{0.0} {}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
	for (auto graph = this->ged_data_.begin(); graph != this->ged_data_.end(); graph++) {
		init_graph_(*graph);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {
	best_feasible_ = NodeMap_(this);
	open_ = std::priority_queue<NodeMap_>();
	omega_ = this->ged_data_.max_edit_cost(g, h) + 10.0;
	Timer timer(time_limit_in_sec_);

	if (not this->initialized_ and lower_bound_method_ == BRANCH_FAST) {
		init_graph_(g);
		init_graph_(h);
	}

	NodeMap_ current_map(g, h, this);
	if (current_map.num_matched_nodes_in_g == g.num_nodes()) {
		extend_half_complete_node_map_(g, h, current_map);
		current_map.matching.set_induced_cost(current_map.induced_cost);
		result.add_node_map(current_map.matching);
		result.set_lower_bound(current_map.induced_cost);
		result.sort_node_maps_and_set_upper_bound();
		return;
	}
	generate_best_child_(g, h, current_map);

	if (not map_root_to_root_) {
		IPFP<UserNodeLabel, UserEdgeLabel> ipfp(this->ged_data_);
		ipfp.set_options("--quadratic-model QAPE --threads " + std::to_string(num_threads_));
		Result ipfp_result;
		ipfp.run_as_util(g, h, ipfp_result);
		best_feasible_.matching = ipfp_result.node_map(0);
		best_feasible_.induced_cost = ipfp_result.upper_bound();
	}

	std::size_t num_itrs{0};
	while (not open_.empty() and not timer.expired()) {
		num_itrs++;
		current_map = open_.top();
		open_.pop();
		if (current_map.lower_bound() >= best_feasible_.induced_cost) {
			continue;
		}
		if ((current_map.num_matched_nodes_in_g == g.num_nodes()) or (current_map.num_matched_nodes_in_h == h.num_nodes())) {
			extend_half_complete_node_map_(g, h, current_map);
			if (current_map.induced_cost < best_feasible_.induced_cost) {
				best_feasible_ = current_map;
			}
			continue;
		}
		if ((current_map.num_matched_nodes_in_g > 0) and current_map.candidates_left()) {
			generate_best_sibling_(g, h, current_map);
		}
		generate_best_child_(g, h, current_map);
	}

	best_feasible_.matching.set_induced_cost(best_feasible_.induced_cost);
	result.add_node_map(best_feasible_.matching);
	result.sort_node_maps_and_set_upper_bound();
	if (open_.empty()) {
		result.set_lower_bound(best_feasible_.induced_cost);
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	lsape_model_ = LSAPESolver::ECBP;
	search_method_ = DFS;
	lower_bound_method_ = BRANCH_FAST;
	num_threads_ = 1;
	time_limit_in_sec_ = 0.0;
	map_root_to_root_ = false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Exact<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	return "[--lsape-model <arg>] [--search-method <arg>] [--lower-bound-method <arg>] [--threads <arg>] [--time-limit] [--map-root-to-root <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Exact<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "lsape-model") {
		if (arg == "EBP") {
			lsape_model_ = LSAPESolver::EBP;
		}
		else if (arg  == "FLWC") {
			lsape_model_ = LSAPESolver::FLWC;
		}
		else if (arg  == "FLCC") {
			lsape_model_ = LSAPESolver::FLCC;
		}
		else if (arg  == "FBP") {
			lsape_model_ = LSAPESolver::FBP;
		}
		else if (arg == "SFBP") {
			lsape_model_ = LSAPESolver::SFBP;
		}
		else if (arg == "FBP0") {
			lsape_model_ = LSAPESolver::FBP0;
		}
		else if (arg  == "ECBP") {
			lsape_model_ = LSAPESolver::ECBP;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option lsape-model. Usage: options = \"[--lsape-model ECBP|EBP|FLWC|FLCC|FBP|SFBP|FBP0] [...]\"");
		}
		return true;
	}
	else if (option == "search-method") {
		if (arg == "DFS") {
			search_method_ = DFS;
		}
		else if (arg  == "ASTAR") {
			search_method_ = ASTAR;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option search-method. Usage: options = \"[--search-method DFS|ASTAR] [...]\"");
		}
		return true;
	}
	else if (option == "lower-bound-method") {
		if (arg == "BRANCH") {
			lower_bound_method_ = BRANCH;
		}
		else if (arg  == "BRANCH_FAST") {
			lower_bound_method_ = BRANCH_FAST;
		}
		else {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option lower-bound-method. Usage: options = \"[--lower-bound-method BRANCH|BRANCH_FAST] [...]\"");
		}
		return true;
	}
	else if (option == "time-limit") {
		try {
			time_limit_in_sec_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option time-limit. Usage: options = \"[--time-limit <convertible to double>] [...]");
		}
		return true;
	}
	else if (option == "threads") {
		try {
			num_threads_ = std::stoi(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument ") + arg + " for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		if (num_threads_ <= 0) {
			throw Error(std::string("Invalid argument ") + arg + " for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		return true;
	}
	else if (option == "map-root-to-root") {
		if (arg == "TRUE") {
			map_root_to_root_ = true;
		}
		else if (arg == "FALSE") {
			map_root_to_root_ = false;
		}
		else {
			throw Error(std::string("Invalid argument ") + arg + " for option map-root-to-root. Usage: options = \"[--map-root-to-root TRUE|FALSE] [...]");
		}
		return true;
	}
	return false;
}

// ==== Definition of private struct Edge_. ====
template<class UserNodeLabel, class UserEdgeLabel>
Exact<UserNodeLabel, UserEdgeLabel>::
Edge_ ::
Edge_(LabelID label, GEDGraph::EdgeID edge_id) :
label{label},
edge_id{edge_id}{}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Exact<UserNodeLabel, UserEdgeLabel>::
Edge_ ::
operator<(const Edge_ & rhs) const {
	return label < rhs.label;
}

// ==== Definition of private class SortedEdges_. ====
template<class UserNodeLabel, class UserEdgeLabel>
Exact<UserNodeLabel, UserEdgeLabel>::
SortedEdges_ ::
SortedEdges_() :
sorted_edges_() {}

template<class UserNodeLabel, class UserEdgeLabel>
Exact<UserNodeLabel, UserEdgeLabel>::
SortedEdges_ ::
SortedEdges_(const GEDGraph & g):
sorted_edges_(){
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		sorted_edges_[*node] = std::vector<Edge_>();
		for (auto edge = g.incident_edges(*node).first; edge != g.incident_edges(*node).second; edge++) {
			sorted_edges_[*node].push_back(Edge_(g.get_edge_label(*edge), *edge));
		}
		std::sort(sorted_edges_[*node].begin(), sorted_edges_[*node].end());
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
SortedEdges_ ::
operator=(const SortedEdges_ & rhs) {
	sorted_edges_ = rhs.sorted_edges_;
}

template<class UserNodeLabel, class UserEdgeLabel>
const std::vector<typename Exact<UserNodeLabel, UserEdgeLabel>::Edge_> &
Exact<UserNodeLabel, UserEdgeLabel>::
SortedEdges_::
get_incident_edges(GEDGraph::NodeID node) const {
	return sorted_edges_.at(node);
}

// ==== Definition of private class NodeMap_. ====
template<class UserNodeLabel, class UserEdgeLabel>
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_ ::
NodeMap_(const GEDGraph & g, const GEDGraph & h, const Exact * exact) :
exact{exact},
matching(),
is_matched_node_in_g(),
is_matched_node_in_h(),
is_candidate_in_h(),
induced_cost{0.0},
lower_bound_to_leaf{0.0},
num_matched_nodes_in_g{0},
num_matched_nodes_in_h{0}{
	for (auto node_g = g.nodes().first; node_g != g.nodes().second; node_g++) {
		is_matched_node_in_g[*node_g] = false;
	}
	for (auto node_h = h.nodes().first; node_h != h.nodes().second; node_h++) {
		is_matched_node_in_h[*node_h] = false;
		is_candidate_in_h[*node_h] = true;
	}
	is_candidate_in_h[GEDGraph::dummy_node()] = true;
	if (exact->map_root_to_root_) {
		num_matched_nodes_in_g++;
		num_matched_nodes_in_h++;
		is_matched_node_in_g[0] = true;
		is_matched_node_in_h[0] = true;
		is_candidate_in_h[0] = false;
		matching.add_assignment(0, 0);
		induced_cost = exact->ged_data_.node_cost(g.get_node_label(0), h.get_node_label(0));
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_ ::
NodeMap_(const NodeMap_ & node_map) :
exact{node_map.exact},
matching(node_map.matching),
is_matched_node_in_g(node_map.is_matched_node_in_g),
is_matched_node_in_h(node_map.is_matched_node_in_h),
is_candidate_in_h(node_map.is_candidate_in_h),
induced_cost{node_map.induced_cost},
lower_bound_to_leaf{node_map.lower_bound_to_leaf},
num_matched_nodes_in_g{node_map.num_matched_nodes_in_g},
num_matched_nodes_in_h{node_map.num_matched_nodes_in_h}{}

template<class UserNodeLabel, class UserEdgeLabel>
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_ ::
NodeMap_(const Exact * exact) :
exact{exact},
matching(),
is_matched_node_in_g(),
is_matched_node_in_h(),
is_candidate_in_h(),
induced_cost{std::numeric_limits<double>::max()},
lower_bound_to_leaf{0.0},
num_matched_nodes_in_g{0},
num_matched_nodes_in_h{0}{}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_ ::
operator=(const NodeMap_ & rhs) {
	matching = rhs.matching;
	is_matched_node_in_g = rhs.is_matched_node_in_g;
	is_matched_node_in_h = rhs.is_matched_node_in_h;
	is_candidate_in_h = rhs.is_candidate_in_h;
	induced_cost = rhs.induced_cost;
	lower_bound_to_leaf = rhs.lower_bound_to_leaf;
	num_matched_nodes_in_g = rhs.num_matched_nodes_in_g;
	num_matched_nodes_in_h = rhs.num_matched_nodes_in_h;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_ ::
lower_bound() const {
	return induced_cost + lower_bound_to_leaf;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_ ::
operator<(const NodeMap_ & rhs) const {
	if (exact->search_method_ == ASTAR){
		return lower_bound() > rhs.lower_bound();
	}
	return num_matched_nodes_in_g < rhs.num_matched_nodes_in_g;
}

template<class UserNodeLabel, class UserEdgeLabel>
GEDGraph::NodeID
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_::
next_unmatched_node_in_g(const GEDGraph & g) const {
	if (num_matched_nodes_in_g < g.num_nodes()) {
		return *(g.nodes().first + num_matched_nodes_in_g);
	}
	return GEDGraph::dummy_node();
}

template<class UserNodeLabel, class UserEdgeLabel>
GEDGraph::NodeID
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_::
last_matched_node_in_g(const GEDGraph & g) const {
	return *(g.nodes().first + num_matched_nodes_in_g - 1);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_::
reset_is_candidate_in_h() {
	//print();
	for (auto is_matched = is_matched_node_in_h.begin(); is_matched != is_matched_node_in_h.end(); is_matched++) {
		is_candidate_in_h[is_matched->first] = not is_matched->second;
	}
	is_candidate_in_h[GEDGraph::dummy_node()] = true;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_::
print() {
	std::cout << "\n\n\t==== node map ====\n\t" << matching;
	std::cout << "\n\t==== matched nodes in h ====\n\t";
	for (auto is_matched = is_matched_node_in_g.begin(); is_matched != is_matched_node_in_g.end(); is_matched++) {
		if (is_matched->second){
			std::cout << is_matched->first << ", ";
		}
	}
	std::cout << "\n\t==== matched nodes in h ====\n\t";
	for (auto is_matched = is_matched_node_in_h.begin(); is_matched != is_matched_node_in_h.end(); is_matched++) {
		if (is_matched->second){
			std::cout << is_matched->first << ", ";
		}
	}
	std::cout << "\n\t==== unmatched nodes in g ====\n\t";
	for (auto is_matched = is_matched_node_in_g.begin(); is_matched != is_matched_node_in_g.end(); is_matched++) {
		if (not is_matched->second){
			std::cout << is_matched->first << ", ";
		}
	}
	std::cout << "\n\t==== unmatched nodes in h ====\n\t";
	for (auto is_matched = is_matched_node_in_h.begin(); is_matched != is_matched_node_in_h.end(); is_matched++) {
		if (not is_matched->second){
			std::cout << is_matched->first << ", ";
		}
	}
	std::cout << "\n\t==== is_candidate_in_h ====\n\t";
	for (auto is_candidate = is_candidate_in_h.begin(); is_candidate != is_candidate_in_h.end(); is_candidate++) {
		if (is_candidate->second){
			std::cout << is_candidate->first << ", ";
		}
	}
	std::cout<<"\n";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Exact<UserNodeLabel, UserEdgeLabel>::
NodeMap_::
candidates_left() {
	for (auto is_candidate = is_candidate_in_h.begin(); is_candidate != is_candidate_in_h.end(); is_candidate++) {
		if (is_candidate->second) {
			return true;
		}
	}
	return false;
}

// ==== Definitions of private helper member functions. ====
template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
init_graph_(const GEDGraph & graph) {
	sorted_edges_[graph.id()] = SortedEdges_(graph);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
append_extension_(const GEDGraph & g, const GEDGraph & h, const NodeMap & extension, NodeMap_ & node_map) {
	std::vector<NodeMap::Assignment> assignments;
	extension.as_relation(assignments);
	for (const auto & assignment : assignments) {
		node_map.matching.add_assignment(assignment.first, assignment.second);
	}
	this->ged_data_.compute_induced_cost(g, h, node_map.matching);
	node_map.induced_cost = node_map.matching.induced_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
init_indices_(const NodeMap_ & node_map, GEDGraph::SizeTNodeMap & g_ids_to_nodes, GEDGraph::SizeTNodeMap & h_ids_to_nodes) const {
	std::size_t id{0};
	for (auto matched_node = node_map.is_matched_node_in_g.begin(); matched_node != node_map.is_matched_node_in_g.end(); matched_node++) {
		if (not matched_node->second) {
			g_ids_to_nodes[id++] = matched_node->first;
		}
	}
	id = 0;
	for (auto matched_node = node_map.is_matched_node_in_h.begin(); matched_node != node_map.is_matched_node_in_h.end(); matched_node++) {
		if (not matched_node->second) {
			h_ids_to_nodes[id++] = matched_node->first;
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
init_master_problem_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & node_map, const GEDGraph::SizeTNodeMap & g_ids_to_nodes, const GEDGraph::SizeTNodeMap & h_ids_to_nodes, DMatrix & master_problem) const {

	// Compute deletion costs.
	for (std::size_t id_i{0}; id_i < master_problem.num_rows() - 1; id_i++) {
		if ((g_ids_to_nodes.at(id_i) == node_map.next_unmatched_node_in_g(g)) and (not node_map.is_candidate_in_h.at(GEDGraph::dummy_node()))) {
			master_problem(id_i, master_problem.num_cols() - 1) = omega_;
		}
		else {
			master_problem(id_i, master_problem.num_cols() - 1) = compute_deletion_cost_(g, node_map, g_ids_to_nodes.at(id_i));
		}
	}

	// Compute insertion costs.
	for (std::size_t id_k{0}; id_k < master_problem.num_cols() - 1; id_k++) {
		master_problem(master_problem.num_rows() - 1, id_k) = compute_insertion_cost_(h, node_map, h_ids_to_nodes.at(id_k));
	}

	// Compute substitution costs in parallel.
#ifdef _OPENMP
	omp_set_num_threads(num_threads_ - 1);
#pragma omp parallel for if(num_threads_ > 1)
#endif
	for (std::size_t row = 0; row < master_problem.num_rows() - 1; row++) {
		for (std::size_t col = 0; col < master_problem.num_cols() - 1; col++) {
			if (lower_bound_method_ == BRANCH) {
				master_problem(row, col) = compute_branch_substitution_cost_(g, h, node_map, g_ids_to_nodes.at(row), h_ids_to_nodes.at(col));
			}
			else {
				master_problem(row, col) = compute_branch_fast_substitution_cost_(g, h, node_map, g_ids_to_nodes.at(row), h_ids_to_nodes.at(col));
			}
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Exact<UserNodeLabel, UserEdgeLabel>::
compute_deletion_cost_(const GEDGraph & g, const NodeMap_ & node_map, GEDGraph::NodeID i) const {
	// Collect node deletion cost.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), ged::dummy_label())};

	// Collect edge deletion costs.
	auto incident_edges_i = g.incident_edges(i);
	for (auto ij = incident_edges_i.first; ij != incident_edges_i.second; ij++) {
		if (not node_map.is_matched_node_in_g.at(g.head(*ij))) {
			cost += this->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label()) * 0.5;
		}
		else {
			cost += this->ged_data_.edge_cost(g.get_edge_label(*ij), ged::dummy_label());
		}
	}

	// Return overall deletion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Exact<UserNodeLabel, UserEdgeLabel>::
compute_insertion_cost_(const GEDGraph & h, const NodeMap_ & node_map, GEDGraph::NodeID k) const {
	// Collect node insertion cost.
	double cost{this->ged_data_.node_cost(ged::dummy_label(), h.get_node_label(k))};

	// Collect edge insertion costs.
	auto incident_edges_k = h.incident_edges(k);
	for (auto kl = incident_edges_k.first; kl != incident_edges_k.second; kl++) {
		if (not node_map.is_candidate_in_h.at(h.head(*kl))) {
			cost += this->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl)) * 0.5;
		}
		else {
			cost += this->ged_data_.edge_cost(ged::dummy_label(), h.get_edge_label(*kl));
		}
	}

	// Return overall insertion cost.
	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
extend_half_complete_node_map_(const GEDGraph & g, const GEDGraph & h, NodeMap_ & node_map) const {
	for (auto matched_node = node_map.is_matched_node_in_g.begin(); matched_node != node_map.is_matched_node_in_g.end(); matched_node++) {
		if (not matched_node->second) {
			node_map.matching.add_assignment(matched_node->first, GEDGraph::dummy_node());
			matched_node->second = true;
			node_map.num_matched_nodes_in_g++;
		}
	}
	for (auto matched_node = node_map.is_matched_node_in_h.begin(); matched_node != node_map.is_matched_node_in_h.end(); matched_node++) {
		if (not matched_node->second) {
			node_map.matching.add_assignment(GEDGraph::dummy_node(), matched_node->first);
			matched_node->second = true;
			node_map.num_matched_nodes_in_h++;
		}
	}
	this->ged_data_.compute_induced_cost(g, h, node_map.matching);
	node_map.induced_cost = node_map.matching.induced_cost();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
generate_next_map_(const GEDGraph & g, const GEDGraph & h, NodeMap_ & next_map, bool update_induced_cost, bool update_upper_bound) {

	// construct LSAPE instance
	GEDGraph::SizeTNodeMap g_ids_to_nodes, h_ids_to_nodes;
	init_indices_(next_map, g_ids_to_nodes, h_ids_to_nodes);
	DMatrix master_problem(g.num_nodes() - next_map.num_matched_nodes_in_g + 1, h.num_nodes() - next_map.num_matched_nodes_in_h + 1, 0.0);
	init_master_problem_(g, h, next_map, g_ids_to_nodes, h_ids_to_nodes, master_problem);

	// solve LSAPE instance and update lower bound to leaf
	LSAPESolver master_problem_solver(master_problem);
	master_problem_solver.set_model(lsape_model_);
	master_problem_solver.solve();
	NodeMap extension;
	util::construct_node_map_from_solver(master_problem_solver, g_ids_to_nodes, h_ids_to_nodes, extension);
	next_map.lower_bound_to_leaf = master_problem_solver.minimal_cost();

	// update matchings
	GEDGraph::NodeID next_node_g{next_map.next_unmatched_node_in_g(g)};
	GEDGraph::NodeID next_node_h{extension.image(next_node_g)};
	next_map.is_matched_node_in_g[next_node_g] = true;
	next_map.matching.add_assignment(next_node_g, next_node_h);
	next_map.num_matched_nodes_in_g++;
	if (next_node_h != GEDGraph::dummy_node()) {
		next_map.is_matched_node_in_h[next_node_h] = true;
		next_map.num_matched_nodes_in_h++;
	}

	// update members is_candidate and induced_cost
	next_map.is_candidate_in_h[next_node_h] = false;
	if (update_induced_cost) {
		update_induced_cost_(g, h, next_map);
	}
	open_.push(next_map);

	if (update_upper_bound) {
		append_extension_(g, h, extension, next_map);
		if (next_map.induced_cost < best_feasible_.induced_cost) {
			best_feasible_ = next_map;
		}
	}

}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
generate_best_child_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & current_map) {
	NodeMap_ child_map(current_map);
	child_map.reset_is_candidate_in_h();
	generate_next_map_(g, h, child_map, true, false);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
generate_best_sibling_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & current_map) {
	NodeMap_ sibling_map(current_map);
	sibling_map.num_matched_nodes_in_g--;
	GEDGraph::NodeID next_node_g{sibling_map.next_unmatched_node_in_g(g)};
	sibling_map.is_matched_node_in_g[next_node_g] = false;
	GEDGraph::NodeID next_node_h{sibling_map.matching.image(next_node_g)};
	if (next_node_h != GEDGraph::dummy_node()) {
		sibling_map.is_matched_node_in_h[next_node_h] = false;
		sibling_map.num_matched_nodes_in_h--;
		sibling_map.matching.erase_pre_image(next_node_h);
	}
	generate_next_map_(g, h, sibling_map, false, false);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Exact<UserNodeLabel, UserEdgeLabel>::
update_induced_cost_(const GEDGraph & g, const GEDGraph & h, NodeMap_ & node_map) const {
	GEDGraph::NodeID i{node_map.last_matched_node_in_g(g)};
	GEDGraph::NodeID k{node_map.matching.image(i)};
	if (k != GEDGraph::dummy_node()) {
		node_map.induced_cost += this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k));
	}
	else {
		node_map.induced_cost += this->ged_data_.node_cost(g.get_node_label(i), dummy_label());
	}
	std::vector<NodeMap::Assignment> assignments;
	node_map.matching.as_relation(assignments);
	for (const auto & assignment : assignments) {
		GEDGraph::NodeID j{assignment.first};
		GEDGraph::NodeID l{assignment.second};
		if (g.is_edge(i, j) and h.is_edge(k, l)) {
			node_map.induced_cost += this->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), h.get_edge_label(h.get_edge(k, l)));
		}
		else if (g.is_edge(i, j)) {
			node_map.induced_cost += this->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), dummy_label());
		}
		else if (h.is_edge(k, l)) {
			node_map.induced_cost += this->ged_data_.edge_cost(dummy_label(), h.get_edge_label(h.get_edge(k, l)));
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Exact<UserNodeLabel, UserEdgeLabel>::
compute_branch_fast_substitution_cost_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & node_map, GEDGraph::NodeID i, GEDGraph::NodeID k) const {
	// Collect node substitution costs.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k))};

	// Collect outer edge costs.
	std::vector<NodeMap::Assignment> assignments;
	node_map.matching.as_relation(assignments);
	for (const auto & assignment : assignments) {
		GEDGraph::NodeID j{assignment.first};
		GEDGraph::NodeID l{assignment.second};
		if (g.is_edge(i, j) and h.is_edge(k, l)) {
			cost += this->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), h.get_edge_label(h.get_edge(k, l)));
		}
		else if (g.is_edge(i, j)) {
			cost += this->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), dummy_label());
		}
		else if (h.is_edge(k, l)) {
			cost += this->ged_data_.edge_cost(dummy_label(), h.get_edge_label(h.get_edge(k, l)));
		}
	}

	// Collect unmatched edge labels.
	std::vector<LabelID> edge_labels_to_unmatched_neighbours_i;
	for (auto ij = sorted_edges_.at(g.id()).get_incident_edges(i).begin(); ij != sorted_edges_.at(g.id()).get_incident_edges(i).end(); ij++) {
		if (not node_map.is_matched_node_in_g.at(g.head(ij->edge_id))) {
			edge_labels_to_unmatched_neighbours_i.push_back(ij->label);
		}
	}
	std::vector<LabelID> edge_labels_to_unmatched_neighbours_k;
	for (auto kl = sorted_edges_.at(h.id()).get_incident_edges(k).begin(); kl != sorted_edges_.at(h.id()).get_incident_edges(k).end(); kl++) {
		if (not node_map.is_matched_node_in_h.at(h.head(kl->edge_id))) {
			edge_labels_to_unmatched_neighbours_k.push_back(kl->label);
		}
	}

	// Compute and add minimal edge insertion costs.
	if (edge_labels_to_unmatched_neighbours_i.size() < edge_labels_to_unmatched_neighbours_k.size()) {
		double min_ins_cost{std::numeric_limits<double>::infinity()};
		for (auto label_h = edge_labels_to_unmatched_neighbours_k.begin(); label_h != edge_labels_to_unmatched_neighbours_k.end(); label_h++) {
			min_ins_cost = std::min(min_ins_cost, this->ged_data_.edge_cost(dummy_label(), *label_h));
		}
		cost += static_cast<double>(edge_labels_to_unmatched_neighbours_k.size() - edge_labels_to_unmatched_neighbours_i.size()) * min_ins_cost * 0.5;
	}

	// Compute and add minimal edge deletion costs.
	if (edge_labels_to_unmatched_neighbours_i.size() > edge_labels_to_unmatched_neighbours_k.size()) {
		double min_del_cost{std::numeric_limits<double>::infinity()};
		for (auto label_g = edge_labels_to_unmatched_neighbours_i.begin(); label_g != edge_labels_to_unmatched_neighbours_i.end(); label_g++) {
			min_del_cost = std::min(min_del_cost, this->ged_data_.edge_cost(*label_g, dummy_label()));
		}
		cost += static_cast<double>(edge_labels_to_unmatched_neighbours_i.size() - edge_labels_to_unmatched_neighbours_k.size()) * min_del_cost * 0.5;
	}

	// Compute minimal edge relabelling costs.
	double min_rel_cost{std::numeric_limits<double>::infinity()};
	for (auto label_g = edge_labels_to_unmatched_neighbours_i.begin(); label_g != edge_labels_to_unmatched_neighbours_i.end(); label_g++) {
		for (auto label_h = edge_labels_to_unmatched_neighbours_k.begin(); label_h != edge_labels_to_unmatched_neighbours_k.end(); label_h++) {
			if (*label_g != *label_h) {
				min_rel_cost = std::min(min_rel_cost, this->ged_data_.edge_cost(*label_g, *label_h));
			}
		}
	}

	// Compute multiset intersection size.
	std::size_t intersection_size{0};
	auto label_g = edge_labels_to_unmatched_neighbours_i.begin();
	auto label_h = edge_labels_to_unmatched_neighbours_k.begin();
	while ((label_g != edge_labels_to_unmatched_neighbours_i.end()) and (label_h != edge_labels_to_unmatched_neighbours_k.end())) {
		if (*label_g == *label_h) {
			intersection_size++;
			label_g++;
			label_h++;
		}
		else if (*label_g < *label_h) {
			label_g++;
		}
		else {
			label_h++;
		}
	}

	// Collect edge relabelling costs.
	std::size_t gamma(std::min(edge_labels_to_unmatched_neighbours_i.size(), edge_labels_to_unmatched_neighbours_k.size()) - intersection_size);
	if (gamma > 0) {
		cost += static_cast<double>(gamma) * min_rel_cost * 0.5;
	}

	return cost;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
Exact<UserNodeLabel, UserEdgeLabel>::
compute_branch_substitution_cost_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & node_map, GEDGraph::NodeID i, GEDGraph::NodeID k) const {
	// Collect node substitution costs.
	double cost{this->ged_data_.node_cost(g.get_node_label(i), h.get_node_label(k))};

	// Collect outer edge costs.
	std::vector<NodeMap::Assignment> assignments;
	node_map.matching.as_relation(assignments);
	for (const auto & assignment : assignments) {
		GEDGraph::NodeID j{assignment.first};
		GEDGraph::NodeID l{assignment.second};
		if (g.is_edge(i, j) and h.is_edge(k, l)) {
			cost += this->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), h.get_edge_label(h.get_edge(k, l)));
		}
		else if (g.is_edge(i, j)) {
			cost += this->ged_data_.edge_cost(g.get_edge_label(g.get_edge(i, j)), dummy_label());
		}
		else if (h.is_edge(k, l)) {
			cost += this->ged_data_.edge_cost(dummy_label(), h.get_edge_label(h.get_edge(k, l)));
		}
	}

	// Initialize subproblem.
	std::vector<LabelID> edge_labels_to_unmatched_neighbours_i;
	for (auto ij = g.incident_edges(i).first; ij != g.incident_edges(i).second; ij++) {
		if (not node_map.is_matched_node_in_g.at(g.head(*ij))) {
			edge_labels_to_unmatched_neighbours_i.push_back(g.get_edge_label(*ij));
		}
	}
	std::vector<LabelID> edge_labels_to_unmatched_neighbours_k;
	for (auto kl = h.incident_edges(k).first; kl != h.incident_edges(k).second; kl++) {
		if (not node_map.is_matched_node_in_h.at(h.head(*kl))) {
			edge_labels_to_unmatched_neighbours_k.push_back(h.get_edge_label(*kl));
		}
	}
	DMatrix subproblem(edge_labels_to_unmatched_neighbours_i.size() + 1, edge_labels_to_unmatched_neighbours_k.size() + 1);

	// Collect edge deletion costs.
	std::size_t row{0};
	for (auto label_ij = edge_labels_to_unmatched_neighbours_i.begin(); label_ij != edge_labels_to_unmatched_neighbours_i.end(); label_ij++, row++) {
		subproblem(row, edge_labels_to_unmatched_neighbours_k.size()) = this->ged_data_.edge_cost(*label_ij, ged::dummy_label()) * 0.5;
	}

	// Collect edge insertion costs.
	std::size_t col{0};
	for (auto label_kl = edge_labels_to_unmatched_neighbours_k.begin(); label_kl != edge_labels_to_unmatched_neighbours_k.end(); label_kl++, col++) {
		subproblem(edge_labels_to_unmatched_neighbours_i.size(), col) = this->ged_data_.edge_cost(ged::dummy_label(), *label_kl) * 0.5;
	}

	// Collect edge relabelling costs.
	row = 0;
	for (auto label_ij = edge_labels_to_unmatched_neighbours_i.begin(); label_ij != edge_labels_to_unmatched_neighbours_i.end(); label_ij++, row++) {
		col = 0;
		for (auto label_kl = edge_labels_to_unmatched_neighbours_k.begin(); label_kl != edge_labels_to_unmatched_neighbours_k.end(); label_kl++, col++) {
			subproblem(row, col) = this->ged_data_.edge_cost(*label_ij, *label_kl) * 0.5;
		}
	}

	// Solve subproblem.
	LSAPESolver subproblem_solver(subproblem);
	subproblem_solver.set_model(lsape_model_);
	subproblem_solver.solve();

	// Update and return overall substitution cost.
	cost += subproblem_solver.minimal_cost();
	return cost;
}

}

#endif /* SRC_METHODS_EXACT_IPP_ */
