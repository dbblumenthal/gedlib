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
 * @file  ged_method.ipp
 * @brief GEDMethod class definition
 */

#ifndef SRC_METHODS_GED_METHOD_IPP_
#define SRC_METHODS_GED_METHOD_IPP_

namespace ged {

// === Definitions of destructors and constructors. ===
template<class UserNodeLabel, class UserEdgeLabel>
GEDMethod<UserNodeLabel, UserEdgeLabel>::
~GEDMethod() {}

template<class UserNodeLabel, class UserEdgeLabel>
GEDMethod<UserNodeLabel, UserEdgeLabel>::
GEDMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
initialized_{false},
ged_data_(ged_data),
options_(),
lower_bound_(0.0),
upper_bound_(std::numeric_limits<double>::infinity()),
node_map_(0, 0),
runtime_(),
init_time_(){}

// Definitions of public member functions.
template<class UserNodeLabel, class UserEdgeLabel>
Seconds
GEDMethod<UserNodeLabel, UserEdgeLabel>::
get_runtime() const {
	return runtime_;
}

template<class UserNodeLabel, class UserEdgeLabel>
Seconds
GEDMethod<UserNodeLabel, UserEdgeLabel>::
get_init_time() const {
	return init_time_;
}

template<class UserNodeLabel, class UserEdgeLabel>
const NodeMap &
GEDMethod<UserNodeLabel, UserEdgeLabel>::
get_node_map() const {
	return node_map_;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDMethod<UserNodeLabel, UserEdgeLabel>::
get_lower_bound() const {
	return lower_bound_;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDMethod<UserNodeLabel, UserEdgeLabel>::
init() {
	auto start = std::chrono::high_resolution_clock::now();
	ged_init_();
	auto end = std::chrono::high_resolution_clock::now();
	init_time_ = end - start;
	initialized_ = true;
}

template<class UserNodeLabel, class UserEdgeLabel>
double
GEDMethod<UserNodeLabel, UserEdgeLabel>::
get_upper_bound() const {
	return upper_bound_;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDMethod<UserNodeLabel, UserEdgeLabel>::
set_options(const std::string & options) {
	util::options_string_to_options_map(options, options_);
	ged_set_default_options_();
	for (auto option_arg : options_) {
		if (not ged_parse_option_(option_arg.first, option_arg.second)) {
			throw Error("Invalid option \"" + option_arg.first + "\". Usage: options = \"" + ged_valid_options_string_() + "\".");
		}
	}
	initialized_ = false;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDMethod<UserNodeLabel, UserEdgeLabel>::
run(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) {
	Result result;
	auto start = std::chrono::high_resolution_clock::now();
	run_as_util(ged_data_.graph(g_id), ged_data_.graph(h_id), result);
	auto end = std::chrono::high_resolution_clock::now();
	lower_bound_ = result.lower_bound();
	upper_bound_ = result.upper_bound();
	if (result.num_node_maps() > 0) {
		node_map_ = result.node_map(0);
	}
	runtime_ = end - start;
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDMethod<UserNodeLabel, UserEdgeLabel>::
run_as_util(const GEDGraph & g, const GEDGraph & h, Result & result) {

	// Compute optimal solution and return if at least one of the two graphs is empty.
	if ((g.num_nodes() == 0) or (h.num_nodes() == 0)) {
		std::size_t index_node_map{result.add_node_map(g.num_nodes(), h.num_nodes())};
		for (auto node = g.nodes().first; node != g.nodes().second; node++) {
			result.node_map(index_node_map).add_assignment(*node, GEDGraph::dummy_node());
		}
		for (auto node = h.nodes().first; node != h.nodes().second; node++) {
			result.node_map(index_node_map).add_assignment(GEDGraph::dummy_node(), *node);
		}
		ged_data_.compute_induced_cost(g, h, result.node_map(index_node_map));
		result.set_lower_bound(result.upper_bound());
		return;
	}

	// Run the method.
	ged_run_(g, h, result);
}

// === Default definitions of private virtual member functions to be overridden by derived classes. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
GEDMethod<UserNodeLabel, UserEdgeLabel>::
ged_init_() {}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDMethod<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {}

template<class UserNodeLabel, class UserEdgeLabel>
bool
GEDMethod<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	return false;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
GEDMethod<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	return "";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
GEDMethod<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {}

}

#endif /* SRC_METHODS_GED_METHOD_IPP_ */
