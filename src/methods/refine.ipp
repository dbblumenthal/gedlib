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
LSBasedMethod<UserNodeLabel, UserEdgeLabel>(ged_data) {}

// === Definitions of member functions inherited from LSBasedMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
Refine<UserNodeLabel, UserEdgeLabel>::
ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) {
	output_node_map = initial_node_map;
	double best_swap_cost{-1.0};
	while ((best_swap_cost < 0) and (output_node_map.induced_cost() > lower_bound)) {
		NodeMap::Assignment best_swap_assignment_1;
		NodeMap::Assignment best_swap_assignment_2;
		best_swap_cost = 0.0;
		for (auto i = g.nodes().first; i != g.nodes().second; i++) {
			for (auto j = g.nodes().first; j != g.nodes().second; j++) {
				if (*i == *j) {
					continue;
				}
				NodeMap::Assignment assignment_1(*i, output_node_map.image(*i));
				NodeMap::Assignment assignment_2(*j, output_node_map.image(*j));
				double swap_cost{this->ged_data_.swap_cost(g, h, assignment_1, assignment_2, output_node_map)};
				if (swap_cost < best_swap_cost) {
					best_swap_cost = swap_cost;
					best_swap_assignment_1 = assignment_1;
					best_swap_assignment_2 = assignment_2;
				}
			}
		}
		if (best_swap_cost < 0) {
			this->ged_data_.swap_assignments(best_swap_assignment_1, best_swap_assignment_2, best_swap_cost, output_node_map);
		}
	}
}

}


#endif /* SRC_METHODS_REFINE_IPP_ */
