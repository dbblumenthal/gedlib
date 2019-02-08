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
 * @file refine.ipp
 * @brief ged::Refine class definition.
 */

#ifndef SRC_METHODS_REFINE_IPP_
#define SRC_METHODS_REFINE_IPP_

namespace ged {

// === Definitions of destructor and constructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
Refine<UserNodeLabel, UserEdgeLabel>::
~Refine() {}

template<class UserNodeLabel, class UserEdgeLabel>
Refine<UserNodeLabel, UserEdgeLabel>::
Refine(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
LSBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
max_swap_size_{2},
naive_{false},
add_dummy_assignment_{true} {}

// === Definitions of member functions inherited from LSBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Refine<UserNodeLabel, UserEdgeLabel>::
ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) {
	output_node_map = initial_node_map;
	double best_swap_cost{0.0};
	Swap_ best_swap;
	std::size_t swap_size{2};
	while ((output_node_map.induced_cost() > lower_bound) and ((best_swap_cost < 0) or (swap_size <= max_swap_size_))) {
		std::vector<NodeMap::Assignment> assignments;
		output_node_map.as_relation(assignments);
		if (add_dummy_assignment_) {
			assignments.emplace_back(GEDGraph::dummy_node(), GEDGraph::dummy_node());
		}
		// terminate if there are not enough assignments to carry out swaps of the current swap size
		if (swap_size > assignments.size()) {
			break;
		}
		// initialization of swapped assignments
		std::vector<std::size_t> swapped_original_indices;
		for (std::size_t i=0; i<swap_size; i++) {
			swapped_original_indices.emplace_back(i);
		}
		// test all possible subsets
		do {
			std::vector<NodeMap::Assignment> swapped_original_assignments;
			for(std::size_t index : swapped_original_indices){
				swapped_original_assignments.emplace_back(assignments[index]);
			}
			// continue if all nodes on the left or all nodes on the right are dummy nodes
			bool real_node_on_left{false};
			bool real_node_on_right{false};
			for (const auto & assignment : swapped_original_assignments) {
				real_node_on_left = (real_node_on_left or (assignment.first != GEDGraph::dummy_node()));
				real_node_on_right = (real_node_on_right or (assignment.second != GEDGraph::dummy_node()));
				if (real_node_on_left and real_node_on_right) {
					break;
				}
			}
			if (not (real_node_on_left and real_node_on_right)) {
				continue;
			}
			// initialization of the swapping cycle inside the current subset
			std::vector<std::size_t> cycle(swap_size-1);
			for (std::size_t i=0; i < swap_size-1;i++) {
				cycle[i] = i+1;
			}
			// test all possible cycle within the swapping set
			do {
				std::vector<NodeMap::Assignment> swapped_new_assignments;
				NodeMap::Assignment original_assignment(swapped_original_assignments.at(0));
				for (std::size_t i=0; i < swap_size - 1; i++) {
					swapped_new_assignments.emplace_back(swapped_original_assignments.at(cycle.at(i)).first, original_assignment.second);
					original_assignment = swapped_original_assignments.at(cycle.at(i));
				}
				swapped_new_assignments.emplace_back(swapped_original_assignments.at(0).first, original_assignment.second);
				Swap_ current_swap;
				current_swap.original_assignments = swapped_original_assignments;
				current_swap.new_assignments = swapped_new_assignments;
				double current_swap_cost = current_swap.cost(g, h, this->ged_data_, output_node_map, naive_);
				if (current_swap_cost - best_swap_cost < -0.000000001) {
					best_swap_cost = current_swap_cost;
					best_swap = current_swap;
				}
			} while (std::next_permutation(cycle.data(), cycle.data() + swap_size-1));
		} while (this->next_subset_(assignments.size(), swapped_original_indices));
		if (best_swap_cost < -0.000000001) {
			best_swap.do_swap(output_node_map, best_swap_cost);
			best_swap_cost = 0.0;
			swap_size = 2;
		}
		else {
			best_swap_cost = 0.0;
			swap_size++;
		}
	}
}


template<class UserNodeLabel, class UserEdgeLabel>
bool
Refine<UserNodeLabel, UserEdgeLabel>::
next_subset_(const std::size_t set_size, std::vector<std::size_t> & subset){
	// Returns FALSE if all subsets have been computed,
	// modifies the in_out subset and returns true otherwise.

	std::size_t subset_size = subset.size();
	std::size_t index_first_changed_element = subset_size - 1;
	while (subset.at(index_first_changed_element) == set_size-subset_size+index_first_changed_element) {
		if (index_first_changed_element==0) {
			return false;
		}
		index_first_changed_element--;
	}
	std::size_t index_first_changed_element_in_set = subset[index_first_changed_element];

	for (std::size_t index_changed_element{index_first_changed_element}; index_changed_element < subset_size; index_changed_element++){
		subset[index_changed_element] = index_first_changed_element_in_set + index_changed_element - index_first_changed_element+1;
	}
	return true;

}

template<class UserNodeLabel, class UserEdgeLabel>
bool
Refine<UserNodeLabel, UserEdgeLabel>::
ls_parse_option_(const std::string & option, const std::string & arg) {
	if (option == "max-swap-size") {
		try {
			max_swap_size_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option max-swap-size. Usage: options = \"[--max-swap-size <convertible to int greater equal 2>]\"");
		}
		if (max_swap_size_ < 2) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option max-swap-size. Usage: options = \"[--max-swap-size <convertible to int greater equal 2>]\"");
		}
		return true;
	}
	else if (option == "naive") {
		if (arg == "TRUE") {
			naive_ = true;
		}
		else if (arg != "FALSE") {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option naive. Usage: options = \"[--naive <TRUE|FALSE>]\"");
		}
		return true;
	}
	else if (option == "add-dummy-assignment") {
		if (arg == "FALSE") {
			add_dummy_assignment_ = false;
		}
		else if (arg != "TRUE") {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option add-dummy-assignment. Usage: options = \"[--add-dummy-assignment <TRUE|FALSE>]\"");
		}
		return true;
	}
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Refine<UserNodeLabel, UserEdgeLabel>::
ls_set_default_options_() {
	max_swap_size_ = 2;
	naive_ = false;
	add_dummy_assignment_ = true;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Refine<UserNodeLabel, UserEdgeLabel>::
ls_valid_options_string_() const {
	return "[--max-swap-size <arg>] [--naive <arg>] [--add-dummy-assignment <arg>]";
}


template<class UserNodeLabel, class UserEdgeLabel>
double
Refine<UserNodeLabel, UserEdgeLabel>::
Swap_::
cost(const GEDGraph & g, const GEDGraph & h, const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data, NodeMap & node_map, bool naive) const {
	double delta{0.0};
	if (naive) {
		double old_induced_cost{node_map.induced_cost()};
		this->do_swap(node_map);
		ged_data.compute_induced_cost(g, h, node_map);
		delta = node_map.induced_cost() - old_induced_cost;
		this->undo_swap(node_map);
		node_map.set_induced_cost(old_induced_cost);
	}
	else if (this->original_assignments.size() == 2){
		delta = ged_data.swap_cost(g, h, original_assignments[0], original_assignments[1], node_map);
	}
	else {
		//collecting swapped vertices in both graphs

		std::vector<GEDGraph::NodeID> g_swapped_vertices;
		std::vector<GEDGraph::NodeID> h_swapped_vertices;

		for(const auto & assignment : this->original_assignments){
			g_swapped_vertices.push_back(assignment.first);
			h_swapped_vertices.push_back(assignment.second);
		}

		//collecting edges incident to swapped vertices in both graphs

		std::vector<GEDGraph::EdgeID> g_incident_edges;
		std::vector<GEDGraph::EdgeID> h_incident_edges;

		for(std::size_t g_vertex_index{0}; g_vertex_index< g_swapped_vertices.size();g_vertex_index++){
			GEDGraph::NodeID g_vertex = g_swapped_vertices[g_vertex_index];
			if (g_vertex != GEDGraph::dummy_node()) {
				for (auto edge = g.incident_edges(g_vertex).first; edge != g.incident_edges(g_vertex).second; edge++) {
					// check if the edge has not been added yet
					bool added_edge = false ;
					for(std::size_t i=0; i<g_vertex_index; ++i){
						if (g.head(*edge) == g_swapped_vertices[i]) {
							added_edge = true;
							break;
						}
					}
					if (!added_edge){
						g_incident_edges.push_back(*edge);
					}
				}
			}
		}
		for(std::size_t h_vertex_index{0}; h_vertex_index< h_swapped_vertices.size();h_vertex_index++){
			GEDGraph::NodeID h_vertex = h_swapped_vertices[h_vertex_index];
			if (h_vertex != GEDGraph::dummy_node()) {
				for (auto edge = h.incident_edges(h_vertex).first; edge != h.incident_edges(h_vertex).second; edge++) {
					// check if the edge has not been added yet
					bool added_edge = false ;
					for(std::size_t i=0; i<h_vertex_index; ++i){
						if (h.head(*edge) == h_swapped_vertices[i]) {
							added_edge = true;
							break;
						}
					}
					if (!added_edge){
						h_incident_edges.push_back(*edge);
					}
				}
			}
		}
		// Compute node cost delta.
		for(auto&& assignment: this->original_assignments){
			delta -= ged_data.node_cost(g.get_node_label(assignment.first), h.get_node_label(assignment.second));
		}
		for(auto&& assignment: this->new_assignments){
			delta += ged_data.node_cost(g.get_node_label(assignment.first), h.get_node_label(assignment.second));
		}
		// Compute negative part of edge cost delta.
		for (const auto & edge : g_incident_edges) {
			delta -= ged_data.edge_cost(g.get_edge_label(edge), h.get_edge_label(node_map.image(g.tail(edge)), node_map.image(g.head(edge))));
		}
		for (const auto & edge : h_incident_edges) {
			if (not g.is_edge(node_map.pre_image(h.tail(edge)), node_map.pre_image(h.head(edge)))) {
				delta -= ged_data.edge_cost(dummy_label(), h.get_edge_label(edge));
			}
		}
		// Carry out the swap, (note : the node_map induced cost is unaffected)
		this->do_swap(node_map);
		// Compute positive part of edge cost delta.
		for (const auto & edge : g_incident_edges) {
			delta += ged_data.edge_cost(g.get_edge_label(edge), h.get_edge_label(node_map.image(g.tail(edge)), node_map.image(g.head(edge))));
		}
		for (const auto & edge : h_incident_edges) {
			if (not g.is_edge(node_map.pre_image(h.tail(edge)), node_map.pre_image(h.head(edge)))) {
				delta += ged_data.edge_cost(dummy_label(), h.get_edge_label(edge));
			}
		}
		// Undo the swap.
		this->undo_swap(node_map);
		// Return the overall swap cost.
	}
	return delta;
}


template<class UserNodeLabel, class UserEdgeLabel>
void
Refine<UserNodeLabel, UserEdgeLabel>::
Swap_::
do_swap(NodeMap & node_map, double delta_cost) const {
	for(auto new_assignment : this->new_assignments){
		node_map.add_assignment(new_assignment.first,new_assignment.second);
	}
	node_map.set_induced_cost(node_map.induced_cost()+delta_cost);
}

template<class UserNodeLabel, class UserEdgeLabel>
void
Refine<UserNodeLabel, UserEdgeLabel>::
Swap_::
undo_swap(NodeMap & node_map) const {
	for(auto new_assignment : this->original_assignments)
		node_map.add_assignment(new_assignment.first,new_assignment.second);
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
Refine<UserNodeLabel, UserEdgeLabel>::
Swap_::
print() const {
	std::stringstream ss;
	ss << "original assignments: { ";
	for (const auto & assignment : original_assignments) {
		ss << assignment << " ";
	}
	ss << "}";
	ss << ", new assignments: { ";
	for (const auto & assignment : new_assignments) {
		ss << assignment << " ";
	}
	ss << "}";
	return ss.str();
}



}


#endif /* SRC_METHODS_REFINE_IPP_ */
