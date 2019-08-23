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
 * @file simulated_annealing.ipp
 * @brief ged::SimulatedAnnealing class definition.
 */

#ifndef SRC_METHODS_SIMULATED_ANNEALING_IPP_
#define SRC_METHODS_SIMULATED_ANNEALING_IPP_

namespace ged {

// === Definitions of constructor and destructor. ===
template<class UserNodeLabel, class UserEdgeLabel>
SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>::
~SimulatedAnnealing() {
	delete lsape_method_;
	delete lower_bound_method_;
}

template<class UserNodeLabel, class UserEdgeLabel>
SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>::
SimulatedAnnealing(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data) :
GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
lsape_method_{new Bipartite<UserNodeLabel, UserEdgeLabel>(this->ged_data_)},
lsape_method_name_("BIPARTITE"),
lsape_method_options_(""),
lower_bound_method_{nullptr},
lower_bound_method_name_("NONE"),
lower_bound_method_options_(""),
num_threads_{1},
num_iterations_{1000},
start_probability_{0.8},
end_probability_{0.01} {}

// === Definitions of member functions inherited from GEDMethod. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {

	DMatrix lsape_instance;
	lsape_method_->populate_instance_and_run_as_util(g, h, result, lsape_instance);

	// Initialize meta parameters.
	double temperature{1.0/std::log(start_probability_)};
	double cooling_factor{1};
	if (num_iterations_ > 1) {
		cooling_factor = std::pow((1.0/std::log(end_probability_)) / temperature, 1.0 / static_cast<double>(num_iterations_ - 1));
	}

	// Initialize lower bound.
	if (lower_bound_method_ and (lower_bound_method_name_ != lsape_method_name_)) {
		Result lower_bound_result;
		lower_bound_method_->run_as_util(g, h, lower_bound_result);
		result.set_lower_bound(std::max(result.lower_bound(), lower_bound_result.lower_bound()));
	}

	// Initialize the order for the candidate generator.
	std::vector<std::size_t> current_order(std::max(g.num_nodes(), h.num_nodes()));
	for (std::size_t pos{0}; pos < current_order.size(); pos++) {
		current_order[pos] = pos;
	}

	// Initialize variables.
	double delta_sum{0};
	NodeMap best_node_map(result.node_maps().at(0));
	NodeMap current_node_map(best_node_map);
	std::size_t num_iterations{0};
	std::size_t num_iterations_without_improvement{0};
	std::random_device rng;
	std::mt19937 urng(rng());
	std::uniform_real_distribution<double> distr(0, 1);

	// Main loop.
	while ((best_node_map.induced_cost() > result.lower_bound()) and (num_iterations++ < num_iterations_)) {
		NodeMap candidate_node_map(g.num_nodes(), h.num_nodes());
		std::vector<std::size_t> candidate_order;
		generate_candidate_(g, h, lsape_instance, current_order, candidate_order, candidate_node_map);
		double delta{std::fabs(candidate_node_map.induced_cost() - current_node_map.induced_cost())};
		delta_sum += delta;
		double delta_avg{delta_sum / static_cast<double>(num_iterations)};
		double random_number{distr(urng)};
		if ((candidate_node_map.induced_cost() < current_node_map.induced_cost()) or (random_number < std::exp((-delta) / (delta_avg * temperature)))) {
			current_node_map = candidate_node_map;
			current_order = candidate_order;
		}
		if (current_node_map.induced_cost() < best_node_map.induced_cost()) {
			best_node_map = current_node_map;
			num_iterations_without_improvement = 0;
		}
		else {
			num_iterations_without_improvement++;
			random_number = distr(urng);
			if (random_number < static_cast<double>(num_iterations_without_improvement) / static_cast<double>(num_iterations_)) {
				std::shuffle(current_order.begin(), current_order.end(), urng);
			}
		}
		temperature *= cooling_factor;
	}

	// Store the best node map in the result if it has lead to an improvement.
	if (not (best_node_map == result.node_maps().at(0))) {
		result.node_maps().emplace(result.node_maps().begin(), best_node_map);
		result.node_maps().pop_back();
	}

}

template<class UserNodeLabel, class UserEdgeLabel>
void
SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
	lsape_method_->init();
	if (lower_bound_method_) {
		lower_bound_method_->init();
	}
}

template<class UserNodeLabel, class UserEdgeLabel>
void
SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
	delete lsape_method_;
	lsape_method_ = new Bipartite<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
	lsape_method_name_ = std::string("BIPARTITE");
	lsape_method_options_ = std::string("");
	delete lower_bound_method_;
	lower_bound_method_ = nullptr;
	lower_bound_method_name_ = std::string("NONE");
	lower_bound_method_options_ = std::string("");
	num_iterations_ = 1000;
	num_threads_ = 1;
	start_probability_ = 0.8;
	end_probability_ = 0.01;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
	return "[--threads <arg>] [--iterations <arg>] [--start-probability <arg>] [--end-probability <arg>] [--lower-bound-method <arg>] [--lsape-method <arg>] [--lsape-method-options <arg>]";
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
	bool is_valid_option{false};
	if (option == "threads") {
		try {
			num_threads_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		if (num_threads_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option threads. Usage: options = \"[--threads <convertible to int greater 0>] [...]");
		}
		if (lsape_method_options_ != "") {
			lsape_method_options_ += " ";
		}
		lsape_method_options_ += "--threads " + std::to_string(num_threads_);
		if (lower_bound_method_options_ != "") {
			lower_bound_method_options_ += " ";
		}
		lower_bound_method_options_ += "--threads " + std::to_string(num_threads_);
		is_valid_option = true;
	}
	else if (option == "iterations") {
		try {
			num_iterations_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option iterations. Usage: options = \"[--iterations <convertible to int greater 0>] [...]");
		}
		if (num_iterations_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg + "\" for option iterations. Usage: options = \"[--iterations <convertible to int greater 0>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "start-probability") {
		try {
			start_probability_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument ") + arg + " for option start-probability. Usage: options = \"[--start-probability <convertible to double between 0 and 1>] [...]");
		}
		if (start_probability_ < 0.0 or start_probability_ > 1.0) {
			throw Error(std::string("Invalid argument ") + arg + " for option start-probability. Usage: options = \"[--start-probability <convertible to double between 0 and 1>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "end-probability") {
		try {
			end_probability_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument ") + arg + " for option end-probability. Usage: options = \"[--end-probability <convertible to double between 0 and 1>] [...]");
		}
		if (end_probability_ < 0.0 or end_probability_ > 1.0) {
			throw Error(std::string("Invalid argument ") + arg + " for option end-probability. Usage: options = \"[--end-probability <convertible to double between 0 and 1>] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "lower-bound-method") {
		lower_bound_method_name_ = arg;
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
	else if (option == "lsape-method") {
		lsape_method_name_ = arg;
		if (arg == "BRANCH_FAST") {
			lsape_method_ = new BranchFast<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BRANCH_UNIFORM") {
			lsape_method_ = new BranchUniform<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BRANCH") {
			lsape_method_ = new Branch<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "NODE") {
			lsape_method_ = new Node<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "RING") {
			lsape_method_ = new Ring<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "SUBGRAPH") {
			lsape_method_ = new Subgraph<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "WALKS") {
			lsape_method_ = new Walks<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "BIPARTITE_ML") {
			lsape_method_ = new BipartiteML<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg == "RINGE_ML") {
			lsape_method_ = new RingML<UserNodeLabel, UserEdgeLabel>(this->ged_data_);
		}
		else if (arg != "BIPARTITE") {
			throw Error("Invalid argument \"" + arg + "\" for option lsape-method. Usage: options = \"[--lsape-method BIPARTITE|BRANCH_FAST|BRANCH_UNIFORM|BRANCH|NODE|RING|SUBGRAPH|WALKS|BIPARTITE_ML|RING_ML] [...]");
		}
		is_valid_option = true;
	}
	else if (option == "lsape-options") {
		std::string lsape_method_options(arg);
		std::size_t bad_option_start{lsape_method_options.find("--threads")};
		std::size_t next_option_start;
		if (bad_option_start != std::string::npos) {
			next_option_start = lsape_method_options.find("--", bad_option_start + 1);
			if (next_option_start != std::string::npos) {
				lsape_method_options = lsape_method_options.substr(0, bad_option_start) + lsape_method_options.substr(next_option_start);
			}
			else {
				lsape_method_options = lsape_method_options.substr(0, bad_option_start);
			}
		}
		if (lsape_method_options_ != "" and lsape_method_options != "") {
			lsape_method_options_ += " ";
		}
		lsape_method_options_ += lsape_method_options;
		is_valid_option = true;
	}
	if (lower_bound_method_) {
		lower_bound_method_->set_options(lower_bound_method_options_);
	}
	lsape_method_->set_options(lsape_method_options_);
	return is_valid_option;
}

// === Definitions of private helper member functions. ===
template<class UserNodeLabel, class UserEdgeLabel>
void
SimulatedAnnealing<UserNodeLabel, UserEdgeLabel>::
generate_candidate_(const GEDGraph & g, const GEDGraph & h, const DMatrix & lsape_instance, const std::vector<std::size_t> & current_order,
		std::vector<std::size_t> & candidate_order, NodeMap & candidate_node_map) const {

	// Slightly change the order.
	candidate_order = current_order;
	std::random_device rng;
	std::mt19937 urng(rng());
	std::uniform_int_distribution<std::size_t> distr(0, candidate_order.size() - 1);
	std::size_t random_pos{distr(urng)};
	std::size_t value_random_pos{candidate_order.at(random_pos)};
	for (std::size_t pos{random_pos}; pos >= 1; pos--) {
		candidate_order[pos] = candidate_order.at(pos - 1);
	}
	candidate_order[0] = value_random_pos;

	// Greedily compute the node map.
	if (g.num_nodes() >= h.num_nodes()) {
		std::vector<bool> is_unassigned_col(h.num_nodes(), true);
		for (std::size_t row : candidate_order) {
			std::size_t best_col{h.num_nodes()};
			for (std::size_t col{0}; col < h.num_nodes(); col++) {
				if (is_unassigned_col.at(col) and (lsape_instance(row, col) < lsape_instance(row, best_col))) {
					best_col = col;
				}
			}
			if (best_col < h.num_nodes()) {
				candidate_node_map.add_assignment(row, best_col);
				is_unassigned_col[best_col] = false;
			}
			else {
				candidate_node_map.add_assignment(row, GEDGraph::dummy_node());
			}
			for (std::size_t col{0}; col < h.num_nodes(); col++) {
				if (is_unassigned_col.at(col)) {
					candidate_node_map.add_assignment(GEDGraph::dummy_node(), col);
				}
			}
		}
	}
	else {
		std::vector<bool> is_unassigned_row(g.num_nodes(), true);
		for (std::size_t col : candidate_order) {
			std::size_t best_row{g.num_nodes()};
			for (std::size_t row{0}; row < g.num_nodes(); row++) {
				if (is_unassigned_row.at(row) and (lsape_instance(row, col) < lsape_instance(best_row, col))) {
					best_row = row;
				}
			}
			if (best_row < g.num_nodes()) {
				candidate_node_map.add_assignment(best_row, col);
				is_unassigned_row[best_row] = false;
			}
			else {
				candidate_node_map.add_assignment(GEDGraph::dummy_node(), col);
			}
			for (std::size_t row{0}; row < g.num_nodes(); row++) {
				if (is_unassigned_row.at(row)) {
					candidate_node_map.add_assignment(row, GEDGraph::dummy_node());
				}
			}
		}
	}
	this->ged_data_.compute_induced_cost(g, h, candidate_node_map);
}

}

#endif /* SRC_METHODS_SIMULATED_ANNEALING_IPP_ */
