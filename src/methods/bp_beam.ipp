/*!
 * @file bp_beam.ipp
 * @brief ged::BPBeam class definition.
 */

#ifndef SRC_METHODS_BP_BEAM_IPP_
#define SRC_METHODS_BP_BEAM_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
~BPBeam() {
	delete ordering_method_;
}

template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
BPBeam(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
beam_size_{5},
num_random_orderings_{0},
ordering_method_{nullptr} {}

// === Definitions of member functions inherited from LSBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BPBeam<UserNodeLabel, UserEdgeLabel>::
ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) {

	// Get the relational representation of the initial node map.
	std::vector<NodeMap::Assignment> assignments;
	initial_node_map.as_relation(assignments);

	// Initialize the output node map.
	output_node_map = initial_node_map;

	// Run BPBeam multiple times if random ordering of the assignments is chosen.
	for (std::size_t iteration{0}; iteration < std::max(std::size_t(1), num_random_orderings_); iteration++) {

		// Order the assignments.
		order_assignments_(g, h, assignments);

		// Initialize the queue open.
		std::vector<TreeNode_> open;
		open.emplace_back(initial_node_map);

		// Iterate as long as open is not empty.
		while (not open.empty()) {

			// Get the tree node top with the cheapest node map and remove it from open.
			TreeNode_ top(open.at(0));
			open.erase(open.begin());

			// Top is a leaf in the search tree, so no children need to be expanded.
			if (top.depth() == assignments.size() - 1) {
				continue;
			}

			// Put all children of top in the queue and update the output node map is a cheaper node map is encountered.
			NodeMap::Assignment assignment_1(assignments.at(top.depth()));
			for (auto assignment_2 = assignments.begin() + top.depth(); assignment_2 != assignments.end(); assignment_2++) {
				open.emplace_back(top, this->ged_data_, g, h, assignment_1, *assignment_2);
				if (open.back().node_map().induced_cost() < output_node_map.induced_cost()) {
					output_node_map = open.back().node_map();
				}
			}

			// Sort the queue and shrink it to the beam size.
			sort_and_shrink_to_beam_size_(open);
		}
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BPBeam<UserNodeLabel, UserEdgeLabel>::
ls_set_default_options_() {
	beam_size_ = 5;
	num_random_orderings_ = 0;
	delete ordering_method_;
	ordering_method_ = nullptr;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BPBeam<UserNodeLabel, UserEdgeLabel>::
ls_parse_option_(const std::string & option, const std::string & arg) {
	bool return_value{false};
	std::string ordering_options("");
	if (option == "beam-size") {
		try {
			beam_size_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option beam-size. Usage: options = \"[--beam-size <convertible to int greater 0>]\"");
		}
		if (beam_size_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option beam-size. Usage: options = \"[--beam-size <convertible to int greater 0>]\"");
		}
		return_value = true;
	}
	else if (option == "ordering-method") {
		if (arg == "RING_ML") {
			ordering_method_ = new RingML<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BIPARTITE_ML") {
			ordering_method_ = new BipartiteML<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "RANDOM") {
			num_random_orderings_ = 5;
		}
		else {
			throw  Error(std::string("Invalid argument \"") + arg + "\" for option ordering-method. Usage: options = \"[--ordering-method RING_ML|BIPARTITE_ML|RANDOM]\"");
		}
		return_value = true;
	}
	else if (option == "ordering-options") {
		ordering_options = arg;
		std::size_t bad_option_start{ordering_options.find("--threads")};
		std::size_t next_option_start;
		if (bad_option_start != std::string::npos) {
			ordering_options = ordering_options.substr(0, bad_option_start);
			next_option_start = ordering_options.find("--", bad_option_start + 1);
			if (next_option_start != std::string::npos) {
				ordering_options += ordering_options.substr(next_option_start);
			}
		}
		bad_option_start = ordering_options.find("--max-num-solutions");
		if (bad_option_start != std::string::npos) {
			ordering_options = ordering_options.substr(0, bad_option_start);
			next_option_start = ordering_options.find("--", bad_option_start + 1);
			if (next_option_start != std::string::npos) {
				ordering_options += ordering_options.substr(next_option_start);
			}
		}
		return_value = true;
	}
	else if (option == "num-random-orderings") {
		try {
			num_random_orderings_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option num-random-orderings. Usage: options = \"[--num-random-orderings <convertible to int greater 0>]\"");
		}
		if (num_random_orderings_ < 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option num-random-orderings. Usage: options = \"[--num-random-orderings <convertible to int greater equal 0>]\"");
		}
		return_value = true;
	}
	if (ordering_method_) {
		if (ordering_options != "") {
			ordering_options += " ";
		}
		ordering_options += "--threads " + std::to_string(this->num_threads_);
		ordering_method_->set_options(ordering_options);
	}
	return return_value;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
BPBeam<UserNodeLabel, UserEdgeLabel>::
ls_valid_options_string_() const {
	return "[--beam-size <arg>] [--ordering-method <arg>] [--ordering-options <arg>] [--num-random-orderings <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BPBeam<UserNodeLabel, UserEdgeLabel>::
ls_init_() {
	if (ordering_method_) {
		ordering_method_->init();
	}
}

// == Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BPBeam<UserNodeLabel, UserEdgeLabel>::
sort_and_shrink_to_beam_size_(std::vector<TreeNode_> & open) const {
	std::sort(open.begin(), open.end());
	while (open.size() > beam_size_) {
		open.pop_back();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
BPBeam<UserNodeLabel, UserEdgeLabel>::
order_assignments_(const GEDGraph & g, const GEDGraph & h, std::vector<NodeMap::Assignment> & assignments) {
	if (ordering_method_) {
		auto compare = [&, this] (const NodeMap::Assignment & lhs, const NodeMap::Assignment & rhs) -> bool {
			double decision_value_lhs{this->ordering_method_->predict(g, h, lhs)};
			double decision_value_rhs{this->ordering_method_->predict(g, h, rhs)};
			return (decision_value_lhs > decision_value_rhs);
		};
		std::sort(assignments.begin(), assignments.end(), compare);
	}
	else if (num_random_orderings_ > 0) {
		std::random_device rng;
		std::mt19937 urng(rng());
		std::shuffle(assignments.begin(), assignments.end(), urng);
	}
}

// == Definitions of private class TreeNode_. ===
template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
TreeNode_(const TreeNode_ & tree_node) :
node_map_(tree_node.node_map_),
depth_{tree_node.depth_} {}

template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
TreeNode_(const NodeMap & node_map) :
node_map_(node_map),
depth_{0} {}

template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
TreeNode_(const TreeNode_ & tree_node, const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data, const GEDGraph & g, const GEDGraph & h, const NodeMap::Assignment & assignment_1, const NodeMap::Assignment & assignment_2) :
node_map_(tree_node.node_map_),
depth_{tree_node.depth_ + 1} {
	double swap_cost{ged_data.swap_cost(g, h, assignment_1, assignment_2, node_map_)};
	ged_data.swap_assignments(assignment_1, assignment_2, swap_cost, node_map_);
}

template<class UserNodeLabel, class UserEdgeLabel>
const NodeMap &
BPBeam<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
node_map() const {
	return node_map_;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::size_t
BPBeam<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
depth() const {
	return depth_;
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
BPBeam<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
operator<(const TreeNode_ & tree_node) const {
	return (node_map_ < tree_node.node_map_);
}

}

#endif /* SRC_METHODS_BP_BEAM_IPP_ */
