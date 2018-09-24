/*!
 * @file ls_based_method.ipp
 * @brief ged::LSBasedMethod class definition.
 */

#ifndef SRC_METHODS_LS_BASED_METHOD_IPP_
#define SRC_METHODS_LS_BASED_METHOD_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
~LSBasedMethod() {
	delete initialization_method_;
	delete lower_bound_method_;
}

template<class UserNodeLabel, class UserEdgeLabel>
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
LSBasedMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
num_threads_{1},
initialization_method_{nullptr},
initialization_options_(""),
lower_bound_method_{nullptr},
lower_bound_method_options_(""),
random_substitution_ratio_{1.0},
num_initial_solutions_{1},
num_runs_from_initial_solutions_{std::numeric_limits<std::size_t>::max()} {}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
	if (initialization_method_) {
		initialization_method_->init();
	}
	if (lower_bound_method_) {
		lower_bound_method_->init();
	}
	ls_init_();
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {

	// Initialize the method for the run between g and h.
	ls_runtime_init_(g, h);

	// Generate the initial node maps and allocate output node maps.
	std::vector<NodeMap> initial_node_maps;
	generate_initial_node_maps(g, h, initial_node_maps, result);
	for (std::size_t node_map_id = 0; node_map_id < initial_node_maps.size(); node_map_id++) {
		result.add_node_map();
	}

	// Initialize lower bound.
	if (lower_bound_method_) {
		Result lower_bound_result;
		lower_bound_method_->run_as_util(g, h, lower_bound_result);
		result.set_lower_bound(std::max(result.lower_bound(), lower_bound_result.lower_bound()));
	}

	// Parallelly run local searches starting at the initial node maps. Stop if the optimal solution has been found or if the desired number of terminated runs has been reached.
	bool found_optimum{false};
	std::size_t terminated_runs{0};
#ifdef _OPENMP
	omp_set_num_threads(num_threads_ - 1);
#pragma omp parallel for if(num_threads_ > 1)
#endif
	for (std::size_t node_map_id = 0; node_map_id < initial_node_maps.size(); node_map_id++) {
		if (not found_optimum and (terminated_runs < num_runs_from_initial_solutions_)) {
			ls_run_from_initial_solution_(g, h, result.lower_bound(), initial_node_maps.at(node_map_id), result.node_map(node_map_id));
#pragma omp critical
			{
				found_optimum = (found_optimum or (result.lower_bound() >= result.node_map(node_map_id).induced_cost()));
				terminated_runs++;
			}
		}
	}

	// Determine the best node map.
	result.sort_node_maps_and_set_upper_bound(num_runs_from_initial_solutions_);

}

template<class UserNodeLabel, class UserEdgeLabel>
bool
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	bool is_valid_option{false};
	if (option == "initialization-method") {
		if (arg == "BIPARTITE_ML") {
			initialization_method_ = new BipartiteML<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BIPARTITE") {
			initialization_method_ = new Bipartite<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BRANCH_FAST") {
			initialization_method_ = new BranchFast<UserNodeLabel, UserEdgeLabel>(this->ged_data_);;
		}
		else if (arg == "BRANCH_UNIFORM") {
			initialization_method_ = new BranchUniform<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BRANCH") {
			initialization_method_ = new Branch<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "NODE") {
			initialization_method_ = new Node<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "RING_ML") {
			initialization_method_ = new RingML<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "RING") {
			initialization_method_ = new Ring<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "SUBGRAPH") {
			initialization_method_ = new Subgraph<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "WALKS") {
			initialization_method_ = new Walks<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg != "RANDOM") {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option initialization-method. Usage: options = \"[--initialization-method BIPARTITE_ML|BIPARTITE|BRANCH_FAST|BRANCH_UNIFORM|BRANCH|NODE|RING_ML|RING|SUBGRAPH|WALKS|RANDOM] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "initialization-options") {
		std::string initialization_options(arg);
		std::size_t bad_option_start{initialization_options.find("--threads")};
		std::size_t next_option_start;
		if (bad_option_start != std::string::npos) {
			next_option_start = initialization_options.find("--", bad_option_start + 1);
			if (next_option_start != std::string::npos) {
				initialization_options = initialization_options.substr(0, bad_option_start) + initialization_options.substr(next_option_start);
			}
			else {
				initialization_options = initialization_options.substr(0, bad_option_start);
			}
		}
		bad_option_start = initialization_options.find("--max-num-solutions");
		if (bad_option_start != std::string::npos) {
			initialization_options = initialization_options.substr(0, bad_option_start);
			next_option_start = initialization_options.find("--", bad_option_start + 1);
			if (next_option_start != std::string::npos) {
				initialization_options = initialization_options.substr(0, bad_option_start) + initialization_options.substr(next_option_start);
			}
			else {
				initialization_options = initialization_options.substr(0, bad_option_start);
			}
		}
		if (initialization_options_ != "") {
			initialization_options_ += " ";
		}
		initialization_options_ += initialization_options;
		is_valid_option = true;
	}
	else if (option == "lower-bound-method") {
		if (arg == "BRANCH") {
			lower_bound_method_ = new Branch<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BRANCH_FAST") {
			lower_bound_method_ = new BranchFast<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BRANCH_TIGHT") {
			lower_bound_method_ = new BranchTight<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
			if (lower_bound_method_options_ != "") {
				lower_bound_method_options_ += " ";
			}
			lower_bound_method_options_ += "--upper-bound NO ";
		}
		else if (arg != "NONE") {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option lower-bound-method. Usage: options = \"[--lower-bound-method BRANCH|BRANCH_FAST|BRANCH_TIGHT|NONE] [...]\"");
		}
		is_valid_option = true;
	}
	else if (option == "random-substitution-ratio") {
		try {
			random_substitution_ratio_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option random-substitution-ratio. Usage: options = \"[--random-substitution-ratio <convertible to double between 0 and 1>]\"");
		}
		if (random_substitution_ratio_ < 0 or random_substitution_ratio_ > 1) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option random-substitution-ratio. Usage: options = \"[--random-substitution-ratio <convertible to double between 0 and 1>]\"");
		}
		is_valid_option = true;
	}
	else if (option == "initial-solutions") {
		try {
			num_initial_solutions_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option num-initial-solutions. Usage: options = \"[--initial-solutions <convertible to int greater 0>]\"");
		}
		if (num_initial_solutions_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option num-initial-solutions. Usage: options = \"[--initial-solutions <convertible to int greater 0>]\"");
		}
		if (initialization_options_ != "") {
			initialization_options_ += " ";
		}
		initialization_options_ += "--max-num-solutions " + std::to_string(num_initial_solutions_);
		is_valid_option = true;
	}
	else if (option == "runs-from-initial-solutions") {
		if (arg == "ALL") {
			num_runs_from_initial_solutions_ = std::numeric_limits<std::size_t>::max();
		}
		else {
			try {
				num_runs_from_initial_solutions_ = std::stoul(arg);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option runs-from-initial-solutions. Usage: options = \"[--runs-from-initial-solutions ALL|<convertible to int greater 0>]\"");
			}
			if (num_runs_from_initial_solutions_ <= 0) {
				throw Error(std::string("Invalid argument \"") + arg + "\" for option runs-from-initial-solutions. Usage: options = \"[--runs-from-initial-solutions ALL|<convertible to int greater 0>]\"");
			}
		}
		is_valid_option = true;
	}
	else if (option == "threads") {
		try {
			num_threads_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>]\"");
		}
		if (num_threads_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>]\"");
		}
		if (initialization_options_ != "") {
			initialization_options_ += " ";
		}
		initialization_options_ += "--threads " + std::to_string(num_threads_);
		if (lower_bound_method_options_ != "") {
			lower_bound_method_options_ += " ";
		}
		lower_bound_method_options_ += "--threads " + std::to_string(num_threads_);
		is_valid_option = true;
	}
	if (initialization_method_) {
		initialization_method_->set_options(initialization_options_);
	}
	if (lower_bound_method_) {
		lower_bound_method_->set_options(lower_bound_method_options_);
	}
	is_valid_option = is_valid_option or ls_parse_option_(option, arg);
	return is_valid_option;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	if (ls_valid_options_string_() == "") {
		return "[--initialization-method <arg>] [--initialization-options <arg>] [--random-substitution-ratio <arg>] [--initial-solutions <arg>] [--runs-from-initial-solutions <arg>] [--threads <arg>]";
	}
	return ls_valid_options_string_() + " [--initialization-method <arg>] [--initialization-options <arg>] [--random-substitution-ratio <arg>] [--initial-solutions <arg>] [--runs-from-initial-solutions <arg>] [--threads <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	delete initialization_method_;
	initialization_method_ = nullptr;
	initialization_options_ = std::string("");
	delete lower_bound_method_;
	lower_bound_method_ = nullptr;
	lower_bound_method_options_ = std::string("");
	random_substitution_ratio_ = 1.0;
	num_initial_solutions_ = 1;
	num_runs_from_initial_solutions_ = std::numeric_limits<std::size_t>::max();
	num_threads_ = 1;
	ls_set_default_options_();
}

// Definitions of private helper member functions.
template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
generate_initial_node_maps(const GEDGraph & g, const GEDGraph & h, std::vector<NodeMap> & initial_node_maps, Result & result) {
	if (initialization_method_) {
		generate_lsape_based_initial_node_maps(g, h, initial_node_maps, result);
	}
	generate_random_initial_node_maps(g, h, initial_node_maps);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
generate_lsape_based_initial_node_maps(const GEDGraph & g, const GEDGraph & h, std::vector<NodeMap> & initial_node_maps, Result & result) {
	Result lsape_result;
	initialization_method_->run_as_util(g, h, lsape_result);
	initial_node_maps = lsape_result.node_maps();
	result.set_lower_bound(lsape_result.lower_bound());

}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
generate_random_initial_node_maps(const GEDGraph & g, const GEDGraph & h, std::vector<NodeMap> & initial_node_maps) {
	std::vector<GEDGraph::NodeID> permutation_g;
	for (auto node = g.nodes().first; node != g.nodes().second; node++) {
		permutation_g.push_back(*node);
	}
	std::vector<GEDGraph::NodeID> permutation_h;
	for (auto node = h.nodes().first; node != h.nodes().second; node++) {
		permutation_h.push_back(*node);
	}
	std::size_t num_substituted_nodes{std::min(g.num_nodes(), h.num_nodes())};
	num_substituted_nodes = std::lround(num_substituted_nodes * random_substitution_ratio_);
	for (std::size_t counter{initial_node_maps.size()}; counter < num_initial_solutions_; counter++) {
		std::random_device rng_g;
		std::mt19937 urng_g(rng_g());
		std::shuffle(permutation_g.begin(), permutation_g.end(), urng_g);
		std::random_device rng_h;
		std::mt19937 urng_h(rng_h());
		std::shuffle(permutation_h.begin(), permutation_h.end(), urng_h);
		initial_node_maps.emplace_back();
		for (auto node = g.nodes().first; node != g.nodes().second; node++) {
			initial_node_maps.back().add_assignment(*node, GEDGraph::dummy_node());
		}
		for (auto node = h.nodes().first; node != h.nodes().second; node++) {
			initial_node_maps.back().add_assignment(GEDGraph::dummy_node(), *node);
		}
		for (std::size_t pos{0}; pos < num_substituted_nodes; pos++) {
			initial_node_maps.back().add_assignment(permutation_g[pos], permutation_h[pos]);
		}
		this->ged_data_.compute_induced_cost(g, h, initial_node_maps.back());
	}
}

// Default definitions of virtual member functions.
template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) {}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ls_init_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ls_runtime_init_(const GEDGraph & g, const GEDGraph & h) {}

template<class UserNodeLabel, class UserEdgeLabel>
bool
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ls_parse_option_(const std::string & option, const std::string & arg) {
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ls_valid_options_string_() const {
	return "";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
LSBasedMethod<UserNodeLabel, UserEdgeLabel>::
ls_set_default_options_() {}

}



#endif /* SRC_METHODS_LS_BASED_METHOD_IPP_ */
