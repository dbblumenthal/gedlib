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
 * @file bp_beam.ipp
 * @brief ged::BPBeam class definition.
 */

#ifndef SRC_METHODS_BP_BEAM_IPP_
#define SRC_METHODS_BP_BEAM_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
~BPBeam() {}

template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
BPBeam(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
beam_size_{5},
num_orderings_{1} {}

// === Definitions of member functions inherited from LSBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
BPBeam<UserNodeLabel, UserEdgeLabel>::
ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) {

	// Get the relational representation of the initial node map.
	std::vector<NodeMap::Assignment> assignments;
	initial_node_map.as_relation(assignments);
	assignments.emplace_back(GEDGraph::dummy_node(), GEDGraph::dummy_node());

	// Initialize the output node map.
	output_node_map = initial_node_map;

	// Run BPBeam multiple times if random ordering of the assignments is chosen.
	for (std::size_t iteration{0}; iteration < num_orderings_; iteration++) {

		// Shuffle the assignments.
		shuffle_assignments_(assignments);

		// Initialize the queue open.
		std::vector<TreeNode_> open;
		open.emplace_back(initial_node_map, assignments);

		// Iterate as long as open is not empty.
		while (not open.empty()) {

			// Get the tree node top with the cheapest node map and remove it from open.
			TreeNode_ top(open.at(0));
			open.erase(open.begin());
			if (top.node_map().induced_cost() < output_node_map.induced_cost()) {
				output_node_map = top.node_map();
			}

			// Top is a leaf in the search tree, so no children need to be expanded.
			if (top.depth() == assignments.size() - 1) {
				continue;
			}

			// Put all children of top in the queue and update the output node map if a cheaper node map is encountered.
			for (std::size_t pos_swapped_assignment{top.depth()}; pos_swapped_assignment < assignments.size(); pos_swapped_assignment++) {
				open.emplace_back(top, this->ged_data_, g, h, pos_swapped_assignment);
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
	num_orderings_ = 1;
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
	else if (option == "num-orderings") {
		try {
			num_orderings_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option num-random-orderings. Usage: options = \"[--num-orderings <convertible to int greater 0>]\"");
		}
		if (num_orderings_ < 1) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option num-random-orderings. Usage: options = \"[--num-orderings <convertible to int greater 0>]\"");
		}
		return_value = true;
	}
	return return_value;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
BPBeam<UserNodeLabel, UserEdgeLabel>::
ls_valid_options_string_() const {
	return "[--beam-size <arg>] [--num-orderings <arg>]";
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
shuffle_assignments_(std::vector<NodeMap::Assignment> & assignments) {
	if (num_orderings_ > 1) {
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
assignments_(tree_node.assignments_),
depth_{tree_node.depth_} {}

template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
TreeNode_(const NodeMap & node_map, const std::vector<NodeMap::Assignment> & assignments) :
node_map_(node_map),
assignments_(assignments),
depth_{0} {}

template<class UserNodeLabel, class UserEdgeLabel>
BPBeam<UserNodeLabel, UserEdgeLabel>::
TreeNode_ ::
TreeNode_(const TreeNode_ & tree_node, const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data, const GEDGraph & g, const GEDGraph & h, std::size_t pos_swapped_assignment) :
node_map_(tree_node.node_map_),
assignments_(tree_node.assignments_),
depth_{tree_node.depth_ + 1} {
	NodeMap::Assignment assignment_1(assignments_.at(depth_ - 1));
	NodeMap::Assignment assignment_2(assignments_.at(pos_swapped_assignment));
	double swap_cost{ged_data.swap_cost(g, h, assignment_1, assignment_2, node_map_)};
	ged_data.swap_assignments(assignment_1, assignment_2, swap_cost, node_map_);
	assignments_[depth_ - 1].second = assignment_2.second;
	assignments_[pos_swapped_assignment].second = assignment_1.second;
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
