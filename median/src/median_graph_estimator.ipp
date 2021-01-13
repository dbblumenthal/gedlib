/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2019 by David B. Blumenthal                              *
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
 * @file median_graph_estimator.ipp
 * @brief @brief ged::MedianGraphEstimator class definition.
 */

#ifndef MEDIAN_SRC_MEDIAN_GRAPH_ESTIMATOR_IPP_
#define MEDIAN_SRC_MEDIAN_GRAPH_ESTIMATOR_IPP_

namespace ged {

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
MedianGraphEstimator(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, bool constant_node_costs):
ged_env_{ged_env},
init_method_{Options::GEDMethod::BRANCH_FAST},
init_options_(""),
descent_method_{Options::GEDMethod::BRANCH_FAST},
descent_options_(""),
refine_method_{Options::GEDMethod::IPFP},
refine_options_(""),
constant_node_costs_{constant_node_costs},
labeled_nodes_{ged_env_->num_node_labels() > 1},
node_del_cost_{ged_env_->node_del_cost(ged_env_->get_node_label(1))},
node_ins_cost_{ged_env_->node_ins_cost(ged_env_->get_node_label(1))},
labeled_edges_{ged_env_->num_edge_labels() > 1},
edge_del_cost_{ged_env_->edge_del_cost(ged_env_->get_edge_label(1))},
edge_ins_cost_{ged_env_->edge_ins_cost(ged_env_->get_edge_label(1))},
init_type_("RANDOM"),
update_order_{true},
num_random_inits_{10},
desired_num_random_inits_{10},
use_real_randomness_{true},
seed_{0},
refine_{true},
time_limit_in_sec_{0},
epsilon_{0.0001},
max_itrs_{100},
max_itrs_without_update_{3},
num_inits_increase_order_{5},
init_type_increase_order_("K-MEANS++"),
max_itrs_increase_order_{10},
print_to_stdout_{2},
median_id_{undefined()},
node_maps_from_median_(),
sum_of_distances_{0},
best_init_sum_of_distances_{std::numeric_limits<double>::infinity()},
converged_sum_of_distances_{std::numeric_limits<double>::infinity()},
runtime_(),
runtime_initialized_(),
runtime_converged_(),
itrs_(),
num_decrease_order_{0},
num_increase_order_{0},
num_converged_descents_{0},
state_{Options::AlgorithmState::TERMINATED} {
	if (ged_env_ == nullptr) {
		throw Error("The environment pointer passed to the constructor of ged::MedianGraphEstimator is null.");
	}
	else if (not ged_env_->initialized()) {
		throw Error("The environment is uninitialized. Call ged::GEDEnv::init() before passing it to the constructor of ged::MedianGraphEstimator.");
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_options(const std::string & options) {
	set_default_options_();
	std::map<std::string, std::string> options_map;
	util::options_string_to_options_map(options, options_map);
	for (const auto & option : options_map) {
		if (option.first == "init-type") {
			init_type_ = option.second;
			if (option.second != "MEDOID" and option.second != "RANDOM" and option.second != "MIN" and option.second != "MAX" and option.second != "MEAN") {
				throw ged::Error(std::string("Invalid argument ") + option.second + " for option init-type. Usage: options = \"[--init-type RANDOM|MEDOID|EMPTY|MIN|MAX|MEAN] [...]\"");
			}
		}
        else if (option.first == "update-order") {
            if (option.second == "TRUE") {
                update_order_ = true;
            }
            else if (option.second == "FALSE") {
                update_order_ = false;
            }
            else {
                throw ged::Error(std::string("Invalid argument ") + option.second + " for option update-order. Usage: options = \"[--update-order TRUE|FALSE] [...]\"");
            }
        }
		else if (option.first == "random-inits") {
			try {
				num_random_inits_ = std::stoul(option.second);
				desired_num_random_inits_ = num_random_inits_;
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option random-inits. Usage: options = \"[--random-inits <convertible to int greater 0>]\"");
			}
			if (num_random_inits_ <= 0) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option random-inits. Usage: options = \"[--random-inits <convertible to int greater 0>]\"");
			}
		}
		else if (option.first == "randomness") {
			if (option.second == "PSEUDO") {
				use_real_randomness_ = false;
			}
			else if (option.second == "REAL") {
				use_real_randomness_ = true;
			}
			else {
				throw Error(std::string("Invalid argument \"") + option.second  + "\" for option randomness. Usage: options = \"[--randomness REAL|PSEUDO] [...]\"");
			}
		}
		else if (option.first == "stdout") {
			if (option.second == "0") {
				print_to_stdout_ = 0;
			}
			else if (option.second == "1") {
				print_to_stdout_ = 1;
			}
			else if (option.second == "2") {
				print_to_stdout_ = 2;
			}
			else {
				throw Error(std::string("Invalid argument \"") + option.second  + "\" for option stdout. Usage: options = \"[--stdout 0|1|2] [...]\"");
			}
		}
		else if (option.first == "refine") {
			if (option.second == "TRUE") {
				refine_ = true;
			}
			else if (option.second == "FALSE") {
				refine_ = false;
			}
			else {
				throw Error(std::string("Invalid argument \"") + option.second  + "\" for option refine. Usage: options = \"[--refine TRUE|FALSE] [...]\"");
			}
		}
		else if (option.first == "time-limit") {
			try {
				time_limit_in_sec_ = std::stod(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option time-limit.  Usage: options = \"[--time-limit <convertible to double>] [...]");
			}
		}
		else if (option.first == "max-itrs") {
			try {
				max_itrs_ = std::stoi(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option max-itrs. Usage: options = \"[--max-itrs <convertible to int>] [...]");
			}
		}
		else if (option.first == "max-itrs-without-update") {
			try {
				max_itrs_without_update_ = std::stoi(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option max-itrs-without-update. Usage: options = \"[--max-itrs-without-update <convertible to int>] [...]");
			}
		}
		else if (option.first == "seed") {
			try {
				seed_ = std::stoul(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option seed. Usage: options = \"[--seed <convertible to int greater equal 0>] [...]");
			}
		}
		else if (option.first == "epsilon") {
			try {
				epsilon_ = std::stod(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option epsilon. Usage: options = \"[--epsilon <convertible to double greater 0>] [...]");
			}
			if (epsilon_ <= 0) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option epsilon. Usage: options = \"[--epsilon <convertible to double greater 0>] [...]");
			}
		}
		else if (option.first == "inits-increase-order") {
			try {
				num_inits_increase_order_ = std::stoul(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option inits-increase-order. Usage: options = \"[--inits-increase-order <convertible to int greater 0>]\"");
			}
			if (num_inits_increase_order_ <= 0) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option inits-increase-order. Usage: options = \"[--inits-increase-order <convertible to int greater 0>]\"");
			}
		}
		else if (option.first == "init-type-increase-order") {
			init_type_increase_order_ = option.second;
			if (option.second != "CLUSTERS" and option.second != "K-MEANS++") {
				throw ged::Error(std::string("Invalid argument ") + option.second + " for option init-type-increase-order. Usage: options = \"[--init-type-increase-order CLUSTERS|K-MEANS++] [...]\"");
			}
		}
		else if (option.first == "max-itrs-increase-order") {
			try {
				max_itrs_increase_order_ = std::stoi(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option max-itrs-increase-order. Usage: options = \"[--max-itrs-increase-order <convertible to int>] [...]");
			}
		}
		else {
			std::string valid_options("[--init-type <arg>] [--update-order <arg>] [--random-inits <arg>] [--randomness <arg>] [--seed <arg>] [--stdout <arg>] ");
			valid_options += "[--time-limit <arg>] [--max-itrs <arg>] [--epsilon <arg>] ";
			valid_options += "[--inits-increase-order <arg>] [--init-type-increase-order <arg>] [--max-itrs-increase-order <arg>]";
			throw Error(std::string("Invalid option \"") + option.first + "\". Usage: options = \"" + valid_options + "\"");
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_init_method(Options::GEDMethod init_method, const std::string & init_options) {
	init_method_ = init_method;
	init_options_ = init_options;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_descent_method(Options::GEDMethod descent_method, const std::string & descent_options) {
	descent_method_ = descent_method;
	descent_options_ = descent_options;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_refine_method(Options::GEDMethod refine_method, const std::string & refine_options) {
	refine_method_ = refine_method;
	refine_options_ = refine_options;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
run(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID median_id) {

	// Sanity checks.
	if (graph_ids.empty()) {
		throw Error("Empty vector of graph IDs, unable to compute median.");
	}
	bool all_graphs_empty{true};
	for (auto graph_id : graph_ids) {
		if (ged_env_->get_num_nodes(graph_id) > 0) {
			all_graphs_empty = false;
			break;
		}
	}
	if (all_graphs_empty) {
		throw Error("All graphs in the collection are empty.");
	}

	// Start timer and record start time.
	auto start = std::chrono::high_resolution_clock::now();
	Timer timer(time_limit_in_sec_);
	median_id_ = median_id;
	state_ = Options::AlgorithmState::TERMINATED;

	// Get ExchangeGraph representations of the input graphs.
	std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> graphs;
	for (auto graph_id : graph_ids) {
		graphs.emplace(graph_id, ged_env_->get_graph(graph_id, true, true, false));
	}

	// Construct initial medians.
	std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> medians;
	construct_initial_medians_(graph_ids, timer, medians);
	auto end_init = std::chrono::high_resolution_clock::now();
	runtime_initialized_ = end_init - start;

	// Reset information about iterations and number of times the median decreases and increases.
	itrs_ = std::vector<std::size_t>(medians.size(), 0);
	num_decrease_order_ = 0;
	num_increase_order_ = 0;
	num_converged_descents_ = 0;

	// Initialize the best median.
	double best_sum_of_distances{std::numeric_limits<double>::infinity()};
	best_init_sum_of_distances_ = std::numeric_limits<double>::infinity();
	std::map<GEDGraph::GraphID, NodeMap> node_maps_from_best_median;
	ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> best_median;

	// Run block gradient descent from all initial medians.
	ged_env_->set_method(descent_method_, descent_options_);
	for (std::size_t median_pos{0}; median_pos < medians.size(); median_pos++) {

		// Terminate if the timer has expired and at least one SOD has been computed.
		if (timer.expired() and (median_pos > 0)) {
			break;
		}

		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			std::cout << "\n===========================================================\n";
			std::cout << "Block gradient descent for initial median " << median_pos + 1 << " of " << medians.size() << ".\n";
			std::cout << "-----------------------------------------------------------\n";
		}

		// Get reference to the median.
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median{medians.at(median_pos)};

		// Load initial median into the environment.
		ged_env_->load_exchange_graph(median, median_id);
		ged_env_->init(ged_env_->get_init_type());

		// Print information about current iteration.
		ged::ProgressBar progress(graph_ids.size());
		if (print_to_stdout_ == 2) {
			std::cout << "\rComputing initial node maps: " << progress << std::flush;
		}

		// Compute node maps and sum of distances for initial median.
		sum_of_distances_ = 0;
		node_maps_from_median_.clear();
		for (auto graph_id : graph_ids) {
			ged_env_->run_method(median_id, graph_id);
			node_maps_from_median_.emplace(graph_id, ged_env_->get_node_map(median_id, graph_id));
			sum_of_distances_ += node_maps_from_median_.at(graph_id).induced_cost();
			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				progress.increment();
				std::cout << "\rComputing initial node maps: " << progress << std::flush;
			}
		}
		best_init_sum_of_distances_ = std::min(best_init_sum_of_distances_, sum_of_distances_);

		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			std::cout << "\n";
		}


		// Run block gradient descent from initial median.
		bool converged{false};
		std::size_t itrs_without_update{0};
		for (; not termination_criterion_met_(converged, timer, itrs_.at(median_pos), itrs_without_update); itrs_[median_pos]++) {

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				std::cout << "\n===========================================================\n";
				std::cout << "Iteration " << itrs_.at(median_pos) + 1 << " for initial median " << median_pos + 1 << " of " << medians.size() << ".\n";
				std::cout << "-----------------------------------------------------------\n";
			}

			// Initialize flags that tell us what happened in the iteration.
			bool median_modified{false};
			bool node_maps_modified{false};
			bool decreased_order{false};
			bool increased_order{false};

			// Update the median.
			median_modified = update_median_(graphs, median);
			if (update_order_ and (not median_modified or itrs_.at(median_pos) == 0)) {
				decreased_order = decrease_order_(graphs, median);
				if (not decreased_order or itrs_.at(median_pos) == 0) {
					increased_order = increase_order_(graphs, median);
				}
			}

			// Update the number of iterations without update of the median.
			if (median_modified or decreased_order or increased_order) {
				itrs_without_update = 0;
			}
			else {
				itrs_without_update++;
			}

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				std::cout << "Loading median to environment: ... " << std::flush;
			}

			// Load the median into the environment.
			ged_env_->load_exchange_graph(median, median_id);

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				std::cout << "done.\n";
				std::cout << "Re-initializing the environment: ... " << std::flush;
			}


			// Re-initialize the environment.
			ged_env_->init(ged_env_->get_init_type());

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				std::cout << "done.\n";
				std::cout << "Updating induced costs: ... " << std::flush;
			}

			// Compute induced costs of the old node maps w.r.t. the updated median.
			for (auto graph_id : graph_ids) {
				ged_env_->compute_induced_cost(median_id, graph_id, node_maps_from_median_.at(graph_id));
			}

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				std::cout << "done.\n";
			}

			// Update the node maps.
			node_maps_modified = update_node_maps_();

			// Update the order of the median if not improvement can be found with the current order.

			// Update the sum of distances.
			double old_sum_of_distances{sum_of_distances_};
			sum_of_distances_ = 0;
			for (const auto & key_val : node_maps_from_median_) {
				sum_of_distances_ += key_val.second.induced_cost();
			}

			// Print information about current iteration.
			if (print_to_stdout_ == 2) {
				std::cout << "Old local SOD: " << old_sum_of_distances << "\n";
				std::cout << "New local SOD: " << sum_of_distances_ << "\n";
				std::cout << "Best converged SOD: " << best_sum_of_distances << "\n";
				std::cout << "Modified median: " << median_modified << "\n";
				std::cout << "Modified node maps: " << node_maps_modified << "\n";
				std::cout << "Decreased order: " << decreased_order << "\n";
				std::cout << "Increased order: " << increased_order << "\n";
				std::cout << "===========================================================\n";
			}

			converged = not (median_modified or node_maps_modified or decreased_order or increased_order);
		}

		// Update the best median.
		if (sum_of_distances_ < best_sum_of_distances) {
			best_sum_of_distances = sum_of_distances_;
			node_maps_from_best_median = node_maps_from_median_;
			best_median = median;
		}

		// Update the number of converged descents.
		if (converged) {
			num_converged_descents_++;
		}
	}

	// Store the best encountered median.
	sum_of_distances_ = best_sum_of_distances;
	node_maps_from_median_ = node_maps_from_best_median;
	ged_env_->load_exchange_graph(best_median, median_id);
	ged_env_->init(ged_env_->get_init_type());
	auto end_descent = std::chrono::high_resolution_clock::now();
	runtime_converged_ = end_descent - start;

	// Refine the sum of distances and the node maps for the converged median.
	converged_sum_of_distances_ = sum_of_distances_;
	if (refine_) {
		improve_sum_of_distances_(timer);
	}

	// Record end time, set runtime and reset the number of initial medians.
	auto end = std::chrono::high_resolution_clock::now();
	runtime_ = end - start;
	num_random_inits_ = desired_num_random_inits_;

	// Print global information.
	if (print_to_stdout_ != 0) {
		std::cout << "\n===========================================================\n";
		std::cout << "Finished computation of generalized median graph.\n";
		std::cout << "-----------------------------------------------------------\n";
		std::cout << "Best SOD after initialization: " << best_init_sum_of_distances_ << "\n";
		std::cout << "Converged SOD: " << converged_sum_of_distances_ << "\n";
		if (refine_) {
			std::cout << "Refined SOD: " << sum_of_distances_ << "\n";
		}
		std::cout << "Overall runtime: " << runtime_.count() << "\n";
		std::cout << "Runtime of initialization: " << runtime_initialized_.count() << "\n";
		std::cout << "Runtime of block gradient descent: " << runtime_converged_.count() - runtime_initialized_.count() << "\n";
		if (refine_) {
			std::cout << "Runtime of refinement: " << runtime_.count() - runtime_converged_.count() << "\n";
		}
		std::cout << "Number of initial medians: " << medians.size() << "\n";
		std::size_t total_itr{0};
		std::size_t num_started_descents{0};
		for (auto itr : itrs_) {
			total_itr += itr;
			if (itr > 0) {
				num_started_descents++;
			}
		}
		std::cout << "Size of graph collection: " << graph_ids.size() << "\n";
		std::cout << "Number of started descents: " << num_started_descents << "\n";
		std::cout << "Number of converged descents: " << num_converged_descents_ << "\n";
		std::cout << "Overall number of iterations: " << total_itr << "\n";
		std::cout << "Overall number of times the order decreased: " << num_decrease_order_<< "\n";
		std::cout << "Overall number of times the order increased: " << num_increase_order_ << "\n";
		std::cout << "===========================================================\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
improve_sum_of_distances_(const Timer & timer) {
	// Use method selected for refinement phase.
	ged_env_->set_method(refine_method_, refine_options_);

	// Print information about current iteration.
	ged::ProgressBar progress(node_maps_from_median_.size());
	if (print_to_stdout_ == 2) {
		std::cout << "\n===========================================================\n";
		std::cout << "\rImproving node maps and SOD for converged median.\n";
		std::cout << "-----------------------------------------------------------\n";
		std::cout << "\rImproving node maps: " << progress << std::flush;
	}

	// Improving the node maps.
	for (auto & key_val : node_maps_from_median_) {
		if (timer.expired()) {
			if (state_ == Options::AlgorithmState::TERMINATED) {
				state_ = Options::AlgorithmState::CONVERGED;
			}
			break;
		}
		GEDGraph::GraphID graph_id{key_val.first};
		NodeMap & node_map{key_val.second};
		ged_env_->run_method(median_id_, graph_id);
		if (ged_env_->get_upper_bound(median_id_, graph_id) < node_map.induced_cost()) {
			node_map = ged_env_->get_node_map(median_id_, graph_id);
		}
		sum_of_distances_ += node_map.induced_cost();
		// Print information.
		if (print_to_stdout_ == 2) {
			progress.increment();
			std::cout << "\rImproving node maps: " << progress << std::flush;
		}
	}
	sum_of_distances_ = 0;
	for (auto & key_val : node_maps_from_median_) {
		sum_of_distances_ += key_val.second.induced_cost();
	}

	// Print information.
	if (print_to_stdout_ == 2) {
		std::cout << "\n===========================================================\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
median_available_() const {
	return (median_id_ != undefined());
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
Options::AlgorithmState
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_state() const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_state().");
	}
	return state_;
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_sum_of_distances(Options::AlgorithmState state) const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_sum_of_distances().");
	}
	if (state == Options::AlgorithmState::INITIALIZED) {
		return best_init_sum_of_distances_;
	}
	if (state == Options::AlgorithmState::CONVERGED) {
		return converged_sum_of_distances_;
	}
	return sum_of_distances_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_distance_from_median(GEDGraph::GraphID graph_id) const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_distance_from_median().");
	}
	if (node_maps_from_median_.find(graph_id) == node_maps_from_median_.end()) {
		throw Error("No distance available for graph with ID " + std::to_string(graph_id) + ".");
	}
	return node_maps_from_median_.at(graph_id).induced_cost();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const NodeMap &
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_node_map_from_median(GEDGraph::GraphID graph_id) const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_node_map_from_median().");
	}
	if (node_maps_from_median_.find(graph_id) == node_maps_from_median_.end()) {
		throw Error("No node map available for graph with ID " + std::to_string(graph_id) + ".");
	}
	return node_maps_from_median_.at(graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const NodeMap &
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_node_map_from_median(GEDGraph::GraphID graph_id) const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling compute_node_map_from_median().");
	}
	ged_env_->run_method(median_id_, graph_id);
	return ged_env_->get_node_map(median_id_, graph_id);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_runtime(Options::AlgorithmState state) const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_runtime().");
	}
	if (state == Options::AlgorithmState::INITIALIZED) {
		return runtime_initialized_.count();
	}
	if (state == Options::AlgorithmState::CONVERGED) {
		return runtime_converged_.count();
	}
	return runtime_.count();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const std::vector<std::size_t> &
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_itrs() const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_num_itrs().");
	}
	return itrs_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_times_order_decreased() const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_num_times_order_decreased().");
	}
	return num_decrease_order_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_times_order_increased() const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_num_times_order_increased().");
	}
	return num_increase_order_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_converged_descents() const {
	if (not median_available_()) {
		throw Error("No median has been computed. Call run() before calling get_num_converged_descents().");
	}
	return num_converged_descents_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> *
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_ged_env() {
	return ged_env_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_default_options_() {
	init_type_ = "RANDOM";
	update_order_ = true;
	num_random_inits_ = 10;
	desired_num_random_inits_ = 10;
	use_real_randomness_ = true;
	seed_ = 0;
	refine_ = true;
	time_limit_in_sec_ = 0;
	epsilon_ = 0.0001;
	max_itrs_ = 100;
	max_itrs_without_update_ = 3;
	num_inits_increase_order_ = 5;
	init_type_increase_order_ = "K-MEANS++";
	max_itrs_increase_order_ = 10;
	print_to_stdout_ = 2;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
construct_initial_medians_(const std::vector<GEDGraph::GraphID> & graph_ids, const Timer & timer, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) {
	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n===========================================================\n";
		std::cout << "Constructing initial median(s).\n";
		std::cout << "-----------------------------------------------------------\n";
	}

	// Compute or sample the initial median(s).
	initial_medians.clear();
	if (init_type_ == "MEDOID") {
		compute_medoid_(graph_ids, timer, initial_medians);
	}
	else if (init_type_ == "MAX") {
		compute_max_order_graph_(graph_ids, initial_medians);
	}
	else if (init_type_ == "MIN") {
		compute_min_order_graph_(graph_ids, initial_medians);
	}
	else if (init_type_ == "MEAN") {
		compute_mean_order_graph_(graph_ids, initial_medians);
	}
	else {
		sample_initial_medians_(graph_ids, initial_medians);
	}

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n===========================================================\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_max_order_graph_(const std::vector<GEDGraph::GraphID> & graph_ids, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) const  {
	GEDGraph::GraphID max_order_id{0};
	std::size_t max_order{0};
	for (auto g_id : graph_ids) {
		if (ged_env_->get_num_nodes(g_id) > max_order) {
			max_order = ged_env_->get_num_nodes(g_id);
			max_order_id = g_id;
		}
	}
	initial_medians.emplace_back(ged_env_->get_graph(max_order_id, true, true, false));
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_min_order_graph_(const std::vector<GEDGraph::GraphID> & graph_ids, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) const  {
	GEDGraph::GraphID min_order_id{0};
	std::size_t min_order{std::numeric_limits<std::size_t>::infinity()};
	for (auto g_id : graph_ids) {
		if (ged_env_->get_num_nodes(g_id) < min_order) {
			min_order = ged_env_->get_num_nodes(g_id);
			min_order_id = g_id;
		}
	}
	initial_medians.emplace_back(ged_env_->get_graph(min_order_id, true, true, false));
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_mean_order_graph_(const std::vector<GEDGraph::GraphID> & graph_ids, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) const  {
	std::vector<std::pair<std::size_t, GEDGraph::GraphID>> order_id_pairs;
	std::size_t sum_orders{0};
	for (auto g_id : graph_ids) {
		sum_orders += ged_env_->get_num_nodes(g_id);
		order_id_pairs.emplace_back(ged_env_->get_num_nodes(g_id), g_id);
	}
	std::sort(order_id_pairs.begin(), order_id_pairs.end());
	double mean_order{static_cast<double>(sum_orders) / static_cast<double>(graph_ids.size())};
	std::size_t median_pos{0};
	for (std::size_t pos{0}; pos < graph_ids.size(); pos++) {
		if (static_cast<double>(order_id_pairs.at(pos).first) >= mean_order) {
			median_pos = pos;
			break;
		}
	}
	if (median_pos > 0) {
		if (mean_order - static_cast<double>(order_id_pairs.at(median_pos - 1).first) < static_cast<double>(order_id_pairs.at(median_pos).first) - mean_order) {
			median_pos--;
		}
	}
	initial_medians.emplace_back(ged_env_->get_graph(order_id_pairs.at(median_pos).second, true, true, false));
}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_medoid_(const std::vector<GEDGraph::GraphID> & graph_ids, const Timer & timer, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) {
	// Use method selected for initialization phase.
	ged_env_->set_method(init_method_, init_options_);

	// Print information about current iteration.
	ged::ProgressBar progress(graph_ids.size());
	if (print_to_stdout_ == 2) {
		std::cout << "\rComputing medoid: " << progress << std::flush;
	}

	// Compute the medoid.
	GEDGraph::GraphID medoid_id{graph_ids.at(0)};
	double best_sum_of_distances{std::numeric_limits<double>::infinity()};
	for (auto g_id : graph_ids) {
		if (timer.expired()) {
			state_ = Options::AlgorithmState::CALLED;
			break;
		}
		double sum_of_distances{0};
		for (auto h_id : graph_ids) {
			ged_env_->run_method(g_id, h_id);
			sum_of_distances += ged_env_->get_upper_bound(g_id, h_id);
		}
		if (sum_of_distances < best_sum_of_distances) {
			best_sum_of_distances = sum_of_distances;
			medoid_id = g_id;
		}
		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			progress.increment();
			std::cout << "\rComputing medoid: " << progress << std::flush;
		}
	}
	initial_medians.emplace_back(ged_env_->get_graph(medoid_id, true, true, false));

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
sample_initial_medians_(const std::vector<GEDGraph::GraphID> & graph_ids, std::vector<ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & initial_medians) {
	// Print information about current iteration.
	ged::ProgressBar progress(num_random_inits_);
	if (print_to_stdout_ == 2) {
		std::cout << "\rSampling initial medians: " << progress << std::flush;
	}

	// Set the number of initial medians to the size of the collection, if it exceeds this size (undone at the end of run()).
	if (num_random_inits_ > graph_ids.size()) {
		num_random_inits_ = graph_ids.size();
	}

	// Set the seed of the random number generator.
	std::mt19937 urng;
	if (use_real_randomness_) {
		std::random_device rng;
		urng.seed(rng());
	}
	else {
		urng.seed(seed_);
	}

	// Sample initial medians.
	std::vector<GEDGraph::GraphID> shuffled_graph_ids(graph_ids);
	for (std::size_t pos{0}; pos < num_random_inits_; pos++) {
		std::shuffle(shuffled_graph_ids.begin(), shuffled_graph_ids.end(), urng);
		initial_medians.emplace_back(ged_env_->get_graph(shuffled_graph_ids.at(pos), true, true, false));
		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			progress.increment();
			std::cout << "\rSampling initial medians: " << progress << std::flush;
		}
	}

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
termination_criterion_met_(bool converged, const Timer & timer, std::size_t itr, std::size_t itrs_without_update) {
	if (timer.expired() or (max_itrs_ >= 0 ? itr >= max_itrs_ : false)) {
		if (state_ == Options::AlgorithmState::TERMINATED) {
			state_ = Options::AlgorithmState::INITIALIZED;
		}
		return true;
	}
	return converged or (max_itrs_without_update_ >= 0 ? itrs_without_update > max_itrs_without_update_ : false);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_median_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) const {

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "Updating median: " << std::flush;
	}

	// Store copy of the old median.
	ged::ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> old_median(median);

	// Update the node labels.
	if (labeled_nodes_) {
		update_node_labels_(graphs, median);
	}

	// Update the edges and their labels.
	update_edges_(graphs, median);

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "done.\n";
	}

	// Return true if the median has changed
	return not (median == old_median);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_node_labels_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) const {

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "nodes ... " << std::flush;
	}

	// Iterate through all nodes of the median.
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		// Collect the labels of the substituted nodes.
		std::vector<UserNodeLabel> node_labels;
		for (const auto & key_val : graphs) {
			GEDGraph::GraphID graph_id{key_val.first};
			const ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & graph{key_val.second};
			GEDGraph::NodeID k{node_maps_from_median_.at(graph_id).image(i)};
			if (k != GEDGraph::dummy_node()) {
				node_labels.emplace_back(graph.node_labels.at(k));
			}
		}

		// Compute the median label and update the median.
		if (node_labels.size() > 0) {
			UserNodeLabel median_label{ged_env_->median_node_label(node_labels)};
			if (ged_env_->node_rel_cost(median.node_labels.at(i), median_label) > epsilon_) {
				median.node_labels[i] = median_label;
			}
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_edges_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) const {

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "edges ... " << std::flush;
	}

	// Clear the adjacency lists of the median and reset number of edges to 0.
	median.num_edges = 0;
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		median.adj_lists.at(i).clear();
	}


	// Iterate through all possible edges (i,j) of the median.
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		for (std::size_t j{i+1}; j < median.num_nodes; j++) {

			// Collect the labels of the edges to which (i,j) is mapped by the node maps.
			std::vector<UserEdgeLabel> edge_labels;
			for (const auto & key_val : graphs) {
				GEDGraph::GraphID graph_id{key_val.first};
				const ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & graph{key_val.second};
				GEDGraph::NodeID k{node_maps_from_median_.at(graph_id).image(i)};
				GEDGraph::NodeID l{node_maps_from_median_.at(graph_id).image(j)};
				if ((k != GEDGraph::dummy_node()) and (l != GEDGraph::dummy_node())) {
					if (graph.adj_matrix.at(k).at(l) == 1) {
						edge_labels.emplace_back(graph.edge_labels.at(std::make_pair(k, l)));
					}
				}
			}

			// Compute the median edge label and the overall edge relabeling cost.
			double rel_cost{0};
			UserEdgeLabel median_label{ged_env_->get_edge_label(1)};
			if (median.adj_matrix.at(i).at(j) == 1) {
				median_label = median.edge_labels.at(std::make_pair(i, j));
			}
			if (labeled_edges_ and (edge_labels.size() > 0)) {
				UserEdgeLabel new_median_label{ged_env_->median_edge_label(edge_labels)};
				if (ged_env_->edge_rel_cost(median_label, new_median_label) > epsilon_) {
					median_label = new_median_label;
				}
				for (const auto & edge_label : edge_labels) {
					rel_cost += ged_env_->edge_rel_cost(median_label, edge_label);
				}
			}

			// Update the median.
			if (rel_cost < (edge_ins_cost_ + edge_del_cost_) * static_cast<double>(edge_labels.size()) - edge_del_cost_ * static_cast<double>(graphs.size())) {
				median.num_edges++;
				median.adj_matrix[i][j] = 1;
				median.adj_matrix[j][i] = 1;
				median.edge_labels[std::make_pair(i, j)] = median_label;
				median.edge_labels[std::make_pair(j, i)] = median_label;
				median.adj_lists.at(i).emplace_back(j, median_label);
				median.adj_lists.at(j).emplace_back(i, median_label);
			}
			else {
				median.adj_matrix[i][j] = 0;
				median.adj_matrix[j][i] = 0;
				if (median.edge_labels.find(std::make_pair(i, j)) != median.edge_labels.end()) {
					median.edge_labels.erase(std::make_pair(i, j));
					median.edge_labels.erase(std::make_pair(j, i));
				}
			}
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_node_maps_() {
	// Print information about current iteration.
	ged::ProgressBar progress(node_maps_from_median_.size());
	if (print_to_stdout_ == 2) {
		std::cout << "\rUpdating node maps: " << progress << std::flush;
	}

	// Update the node maps.
	bool node_maps_were_modified{false};
	for (auto & key_val : node_maps_from_median_) {
		GEDGraph::GraphID graph_id{key_val.first};
		NodeMap & node_map{key_val.second};
		ged_env_->run_method(median_id_, graph_id);
		if (ged_env_->get_upper_bound(median_id_, graph_id) < node_map.induced_cost() - epsilon_) {
			node_map = ged_env_->get_node_map(median_id_, graph_id);
			node_maps_were_modified = true;
		}
		// Print information about current iteration.
		if (print_to_stdout_ == 2) {
			progress.increment();
			std::cout << "\rUpdating node maps: " << progress << std::flush;
		}
	}

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "\n";
	}

	// Return true if the node maps were modified.
	return node_maps_were_modified;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
decrease_order_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) {

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "Trying to decrease order: ... " << std::flush;
	}

	// Initialize ID of the node that is to be deleted.
	std::size_t id_deleted_node{undefined()};
	bool decreased_order{false};

	// Decrease the order as long as the best deletion delta is negative.
	while (compute_best_deletion_delta_(graphs, median, id_deleted_node) < -epsilon_) {
		decreased_order = true;
		delete_node_from_median_(id_deleted_node, median);
	}

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "done.\n";
	}

	// Return true iff the order was decreased.
	return decreased_order;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_best_deletion_delta_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, const ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median, std::size_t & id_deleted_node) const {
	double best_delta{0};

	// Determine node that should be deleted (if any).
	for (std::size_t i{0}; i < median.num_nodes; i++) {

		// Compute cost delta.
		double delta{0};
		for (const auto & key_val : graphs) {
			GEDGraph::GraphID graph_id{key_val.first};
			const ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & graph{key_val.second};
			GEDGraph::NodeID k{node_maps_from_median_.at(graph_id).image(i)};
			if (k == GEDGraph::dummy_node()) {
				delta -= node_del_cost_;
			}
			else {
				delta += node_ins_cost_ - ged_env_->node_rel_cost(median.node_labels.at(i), graph.node_labels.at(k));
			}
			for (const auto & j_label : median.adj_lists.at(i)) {
				GEDGraph::NodeID l{node_maps_from_median_.at(graph_id).image(j_label.first)};
				if ((k == GEDGraph::dummy_node()) or (l == GEDGraph::dummy_node())) {
					delta -= edge_del_cost_;
				}
				else if (graph.adj_matrix.at(k).at(l) == 0) {
					delta -= edge_del_cost_;
				}
				else {
					delta += edge_ins_cost_ - ged_env_->edge_rel_cost(j_label.second, graph.edge_labels.at(std::make_pair(k, l)));
				}
			}
		}

		// Update best deletion delta.
		if (delta < best_delta - epsilon_) {
			best_delta = delta;
			id_deleted_node = i;
		}
	}

	// Return the best deletion delta.
	return best_delta;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
delete_node_from_median_(std::size_t id_deleted_node, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) {
	// Update the nodes of the median.
	median.num_nodes--;
	median.node_labels.erase(median.node_labels.begin() + id_deleted_node);
	median.original_node_ids.erase(median.original_node_ids.begin() + id_deleted_node);

	// Update the edges of the median.
	median.num_edges = 0;
	median.adj_lists.erase(median.adj_lists.begin() + id_deleted_node);
	for (auto & adj_list : median.adj_lists) {
		adj_list.clear();
	}
	std::vector<std::vector<std::size_t>> new_adj_matrix(median.num_nodes, std::vector<std::size_t>(median.num_nodes, 0));
	std::map<std::pair<std::size_t, std::size_t>, UserEdgeLabel> new_edge_labels;
	for (std::size_t i{0}; i <= median.num_nodes; i++) {
		if (i != id_deleted_node) {
			std::size_t new_i{i < id_deleted_node ? i : i - 1};
			for (std::size_t j{i+1}; j < median.num_nodes; j++) {
				if (j != id_deleted_node) {
					std::size_t new_j{j < id_deleted_node ? j : j - 1};
					if (median.adj_matrix.at(i).at(j) == 1) {
						UserEdgeLabel label{median.edge_labels.at(std::make_pair(i,j))};
						median.num_edges++;
						new_adj_matrix[new_i][new_j] = 1;
						new_adj_matrix[new_j][new_i] = 1;
						new_edge_labels[std::make_pair(new_i, new_j)] = label;
						new_edge_labels[std::make_pair(new_j, new_i)] = label;
						median.adj_lists.at(new_i).emplace_back(new_j, label);
						median.adj_lists.at(new_j).emplace_back(new_i, label);
					}
				}
			}
		}
	}
	median.adj_matrix = new_adj_matrix;
	median.edge_labels = new_edge_labels;

	// Update the node maps.
	for (auto & key_val : node_maps_from_median_) {
		NodeMap & node_map{key_val.second};
		NodeMap new_node_map(median.num_nodes, node_map.num_target_nodes());
		std::vector<bool> is_unassigned_target_node(node_map.num_target_nodes(), true);
		for (std::size_t i{0}; i <= median.num_nodes; i++) {
			if (i != id_deleted_node) {
				std::size_t new_i{i < id_deleted_node ? i : i - 1};
				std::size_t k{node_map.image(i)};
				new_node_map.add_assignment(new_i, k);
				if (k != GEDGraph::dummy_node()) {
					is_unassigned_target_node[k] = false;
				}
			}
		}
		for (std::size_t k{0}; k < node_map.num_target_nodes(); k++) {
			if (is_unassigned_target_node.at(k)) {
				new_node_map.add_assignment(GEDGraph::dummy_node(), k);
			}
		}
		node_map = new_node_map;
	}

	// Increase overall number of decreases.
	num_decrease_order_++;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
increase_order_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) {

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "Trying to increase order: ... " << std::flush;
	}

	// Initialize the best configuration and the best label of the node that is to be inserted.
	std::map<GEDGraph::GraphID, std::size_t> best_config;
	UserNodeLabel best_label{ged_env_->get_node_label(1)};
	bool increased_order{false};

	// Increase the order as long as the best insertion delta is negative.
	while (compute_best_insertion_delta_(graphs, best_config, best_label) < -epsilon_) {
		increased_order = true;
		add_node_to_median_(best_config, best_label, median);
	}

	// Print information about current iteration.
	if (print_to_stdout_ == 2) {
		std::cout << "done.\n";
	}

	// Return true iff the order was increased.
	return increased_order;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_best_insertion_delta_(const std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & graphs, 
		std::map<GEDGraph::GraphID, std::size_t> & best_config, UserNodeLabel & best_label) const {

	// Construct sets of inserted nodes.
	bool no_inserted_node{true};
	std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> inserted_nodes;
	for (const auto & key_val : graphs) {
		GEDGraph::GraphID graph_id{key_val.first};
		const ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & graph{key_val.second};
		inserted_nodes.emplace(graph_id, std::vector<std::pair<std::size_t, UserNodeLabel>>());
		best_config.emplace(graph_id, GEDGraph::dummy_node());
		for (std::size_t k{0}; k < graph.num_nodes; k++) {
			if (node_maps_from_median_.at(graph_id).pre_image(k) == GEDGraph::dummy_node()) {
				no_inserted_node = false;
				inserted_nodes.at(graph_id).emplace_back(std::make_pair(k, graph.node_labels.at(k)));
			}
		}
	}

	// Return 0.0 if no node is inserted in any of the graphs.
	if (no_inserted_node) {
		return 0.0;
	}

	// Compute insertion configuration, label, and delta.
	double best_delta;
	if (not labeled_nodes_) {
		best_delta = compute_insertion_delta_unlabeled_(inserted_nodes, best_config, best_label);
	}
	else if (constant_node_costs_) {
		best_delta = compute_insertion_delta_constant_(inserted_nodes, best_config, best_label);
	}
	else {
		best_delta = compute_insertion_delta_generic_(inserted_nodes, best_config, best_label);
	}

	// Return the best delta.
	return best_delta;

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_insertion_delta_unlabeled_(const std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> & inserted_nodes,
		std::map<GEDGraph::GraphID, std::size_t> & best_config, UserNodeLabel & best_label) const {
	// Construct the nest configuration and compute its insertion delta.
	double best_delta{0};
	best_config.clear();
	for (const auto & key_val : inserted_nodes) {
		GEDGraph::GraphID graph_id{key_val.first};
		const std::vector<std::pair<std::size_t, UserNodeLabel>> & node_set{key_val.second};
		if (node_set.empty()) {
			best_config[graph_id] = GEDGraph::dummy_node();
			best_delta += node_del_cost_;
		}
		else {
			best_config[graph_id] = node_set.at(0).first;
			best_delta -= node_ins_cost_;
		}
	}

	// Return the best insertion delta.
	return best_delta;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_insertion_delta_constant_(const std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> & inserted_nodes,
		std::map<GEDGraph::GraphID, std::size_t> & best_config, UserNodeLabel & best_label) const {

	// Construct histogram and inverse label maps.
	std::map<UserNodeLabel, std::size_t> hist;
	std::map<GEDGraph::GraphID, std::map<UserNodeLabel, std::size_t>> inverse_label_maps;
	for (const auto & key_val : inserted_nodes) {
		GEDGraph::GraphID graph_id{key_val.first};
		const std::vector<std::pair<std::size_t, UserNodeLabel>> & node_set{key_val.second};
		inverse_label_maps.emplace(graph_id, std::map<UserNodeLabel, std::size_t>());
		for (const auto & node : node_set) {
			std::size_t k{node.first};
			UserNodeLabel label{node.second};
			if (inverse_label_maps.at(graph_id).find(label) == inverse_label_maps.at(graph_id).end()) {
				inverse_label_maps.at(graph_id)[label] = k;
				if (hist.find(label) == hist.end()) {
					hist[label] = 1;
				}
				else {
					hist[label]++;
				}
			}
		}
	}

	// Determine the best label.
	std::size_t best_count{0};
	for (const auto & key_val : hist) {
		if (key_val.second > best_count) {
			best_count = key_val.second;
			best_label = key_val.first;
		}
	}

	// Construct the best configuration and compute its insertion delta.
	best_config.clear();
	double best_delta{0};
	double node_rel_cost{ged_env_->node_rel_cost(ged_env_->get_node_label(1), ged_env_->get_node_label(2))};
	bool triangle_ineq_holds{node_rel_cost <= node_del_cost_ + node_ins_cost_};
	for (const auto & key_val : inserted_nodes) {
		GEDGraph::GraphID graph_id{key_val.first};
		if (inverse_label_maps.at(graph_id).find(best_label) != inverse_label_maps.at(graph_id).end()) {
			best_config[graph_id] = inverse_label_maps.at(graph_id).at(best_label);
			best_delta -= node_ins_cost_;
		}
		else if (triangle_ineq_holds and (not inserted_nodes.at(graph_id).empty())) {
			best_config[graph_id] = inserted_nodes.at(graph_id).at(0).first;
			best_delta += node_rel_cost - node_ins_cost_;
		}
		else {
			best_config[graph_id] = GEDGraph::dummy_node();
			best_delta += node_del_cost_;
		}
	}

	// Return the best insertion delta.
	return best_delta;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_insertion_delta_generic_(const std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> & inserted_nodes,
		std::map<GEDGraph::GraphID, std::size_t> & best_config, UserNodeLabel & best_label) const {

	// Collect all node labels of inserted nodes.
	std::vector<UserNodeLabel> node_labels;
	for (const auto & key_val : inserted_nodes) {
		const std::vector<std::pair<std::size_t, UserNodeLabel>> & node_set{key_val.second};
		for (const auto & node : node_set) {
			node_labels.emplace_back(node.second);
		}
	}

	// Compute node label medians that serve as initial solutions for block gradient descent.
	std::vector<UserNodeLabel> initial_node_labels;
	compute_initial_node_labels_(node_labels, initial_node_labels);

	// Determine best insertion configuration, label, and delta via parallel block gradient descent from all initial node labels.
	double best_delta{0};
	for (UserNodeLabel & node_label : initial_node_labels) {
		// Construct local configuration.
		std::map<GEDGraph::GraphID, std::pair<std::size_t, UserNodeLabel>> config;
		for (const auto & key_val : inserted_nodes) {
			GEDGraph::GraphID graph_id{key_val.first};
			config.emplace(graph_id, std::make_pair(GEDGraph::dummy_node(), ged_env_->get_node_label(1)));
		}

		// Run block gradient descent.
		bool converged{false};
		for (std::size_t itr{0}; not insertion_termination_criterion_met_(converged, itr); itr++) {
			converged = (not update_config_(node_label, inserted_nodes, config, node_labels));
			converged = converged and (not update_node_label_(node_labels, node_label));
		}

		// Compute insertion delta of converged solution.
		double delta{0};
		for (const auto & key_val : config) {
			const std::pair<std::size_t, UserNodeLabel> & node{key_val.second};
			if (node.first == GEDGraph::dummy_node()) {
				delta += node_del_cost_;
			}
			else {
				delta += ged_env_->node_rel_cost(node_label, node.second) - node_ins_cost_;
			}
		}

		// Update best delta and global configuration if improvement has been found.
		if (delta < best_delta - epsilon_) {
			best_delta = delta;
			best_label = node_label;
			best_config.clear();
			for (const auto & key_val : config) {
				GEDGraph::GraphID graph_id{key_val.first};
				std::size_t k{key_val.second.first};
				best_config[graph_id] = k;
			}
		}
	}

	// Return the best delta.
	return best_delta;

}


template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_initial_node_labels_(const std::vector<UserNodeLabel> & node_labels, std::vector<UserNodeLabel> & median_labels) const {

	median_labels.clear();
	std::mt19937 urng;
	if (use_real_randomness_) {
		std::random_device rng;
		urng.seed(rng());
	}
	else {
		urng.seed(seed_);
	}

	// Generate the initial node label medians.
	if (init_type_increase_order_ == "K-MEANS++") {

		// Use k-means++ heuristic to generate the initial node label medians.
		std::vector<bool> already_selected(node_labels.size(), false);
		std::uniform_int_distribution<std::size_t> unif_dist(0, node_labels.size() - 1);
		std::size_t selected_label_id{unif_dist(urng)};
		median_labels.emplace_back(node_labels.at(selected_label_id));
		already_selected[selected_label_id] = true;
		while (median_labels.size()  < num_inits_increase_order_) {
			std::vector<double> weights(node_labels.size(), std::numeric_limits<double>::infinity());
			for (std::size_t label_id{0}; label_id < node_labels.size(); label_id++) {
				if (already_selected.at(label_id)) {
					weights[label_id] = 0;
					continue;
				}
				for (const auto & label : median_labels) {
					weights[label_id] = std::min(weights.at(label_id), ged_env_->node_rel_cost(label, node_labels.at(label_id)));
				}
			}
			std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
			selected_label_id = dist(urng);
			median_labels.emplace_back(node_labels.at(selected_label_id));
			already_selected[selected_label_id] = true;
		}
	}
	else {

		// Compute the initial node medians as the medians of randomly generated clusters of (roughly) equal size.
		std::vector<UserNodeLabel> shuffled_node_labels(node_labels);
		std::shuffle(shuffled_node_labels.begin(), shuffled_node_labels.end(), urng);
		std::size_t cluster_size{node_labels.size() / num_inits_increase_order_};
		std::size_t pos{0};
		std::vector<UserNodeLabel> cluster;
		while (median_labels.size()  < num_inits_increase_order_ - 1) {
			while (pos++ < (median_labels.size() + 1) * cluster_size) {
				cluster.emplace_back(shuffled_node_labels.at(pos));
			}
			median_labels.emplace_back(ged_env_->median_node_label(cluster));
			cluster.clear();
		}
		while (pos++ < shuffled_node_labels.size()) {
			cluster.emplace_back(shuffled_node_labels.at(pos));
		}
		median_labels.emplace_back(ged_env_->median_node_label(cluster));
		cluster.clear();
	}

	// Run Lloyd's Algorithm.
	bool converged{false};
	std::vector<std::size_t> closest_median_ids(node_labels.size(), undefined());
	std::vector<std::vector<UserNodeLabel>> clusters(median_labels.size());
	for (std::size_t itr{1}; not insertion_termination_criterion_met_(converged, itr); itr++) {
		converged = not update_clusters_(node_labels, median_labels, closest_median_ids);
		if (not converged) {
			for (auto & cluster : clusters) {
				cluster.clear();
			}
			for (std::size_t label_id{0}; label_id < node_labels.size(); label_id++) {
				clusters.at(closest_median_ids.at(label_id)).emplace_back(node_labels.at(label_id));
			}
			for (std::size_t cluster_id{0}; cluster_id < clusters.size(); cluster_id++) {
				update_node_label_(clusters.at(cluster_id), median_labels.at(cluster_id));
			}
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
insertion_termination_criterion_met_(bool converged, std::size_t itr) const {
	return converged or (max_itrs_increase_order_ >= 0 ? itr >= max_itrs_increase_order_ : false);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_config_(const UserNodeLabel & node_label, const std::map<GEDGraph::GraphID, std::vector<std::pair<std::size_t, UserNodeLabel>>> & inserted_nodes, std::map<GEDGraph::GraphID, std::pair<std::size_t, UserNodeLabel>> & config, std::vector<UserNodeLabel> & node_labels) const {

	// Determine the best configuration.
	bool config_modified{false};
	for (const auto & key_val : inserted_nodes) {
		GEDGraph::GraphID graph_id{key_val.first};
		std::pair<std::size_t, UserNodeLabel> best_assignment{config.at(graph_id)};
		double best_cost{0};
		if (best_assignment.first == GEDGraph::dummy_node()) {
			best_cost = node_del_cost_;
		}
		else {
			best_cost = ged_env_->node_rel_cost(node_label, best_assignment.second) - node_ins_cost_;
		}
		const std::vector<std::pair<std::size_t, UserNodeLabel>> & node_set{key_val.second};
		for (const auto & node : node_set) {
			double cost{ged_env_->node_rel_cost(node_label, node.second) - node_ins_cost_};
			if (cost < best_cost - epsilon_) {
				best_cost = cost;
				best_assignment = node;
				config_modified = true;
			}
		}
		if (node_del_cost_ < best_cost - epsilon_) {
			best_cost = node_del_cost_;
			best_assignment.first = GEDGraph::dummy_node();
			config_modified = true;
		}
		config[graph_id] = best_assignment;
	}

	// Collect the node labels contained in the best configuration.
	node_labels.clear();
	for (const auto & key_val : config) {
		if (key_val.second.first != GEDGraph::dummy_node()) {
			node_labels.emplace_back(key_val.second.second);
		}
	}

	// Return true if the configuration was modified.
	return config_modified;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_node_label_(const std::vector<UserNodeLabel> & node_labels, UserNodeLabel & node_label) const {
	UserNodeLabel new_node_label{ged_env_->median_node_label(node_labels)};
	if (ged_env_->node_rel_cost(new_node_label, node_label) > epsilon_) {
		node_label = new_node_label;
		return true;
	}
	return false;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_clusters_(const std::vector<UserNodeLabel> & node_labels, const std::vector<UserNodeLabel> & median_labels, std::vector<std::size_t> & closest_median_ids) const {

	// Determine the closest median for each node label.
	bool clusters_modified{false};
	for (std::size_t label_id{0}; label_id < node_labels.size(); label_id++) {
		std::size_t closest_median_id{undefined()};
		double dist_to_closest_median{std::numeric_limits<double>::infinity()};
		for (std::size_t median_id{0}; median_id < median_labels.size(); median_id++) {
			double dist_to_median{ged_env_->node_rel_cost(median_labels.at(median_id), node_labels.at(label_id))};
			if (dist_to_median < dist_to_closest_median - epsilon_) {
				dist_to_closest_median = dist_to_median;
				closest_median_id = median_id;
			}
		}
		if (closest_median_id != closest_median_ids.at(label_id)) {
			closest_median_ids[label_id] = closest_median_id;
			clusters_modified = true;
		}
	}

	// Return true if the clusters were modified.
	return clusters_modified;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel>::
add_node_to_median_(const std::map<GEDGraph::GraphID, std::size_t> & best_config, const UserNodeLabel & best_label,
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median) {
	// Update the median.
	for (std::size_t i{0}; i < median.num_nodes; i++) {
		median.adj_matrix.at(i).emplace_back(0);
	}
	median.num_nodes++;
	median.node_labels.emplace_back(best_label);
	median.original_node_ids.emplace_back(util::new_string_or_numeric_(median.original_node_ids));
	median.adj_matrix.emplace_back(std::vector<std::size_t>(median.num_nodes, 0));
	median.adj_lists.emplace_back(std::list<std::pair<std::size_t, UserEdgeLabel>>());

	// Update the node maps.
	for (auto & key_val : node_maps_from_median_) {
		GEDGraph::GraphID graph_id{key_val.first};
		NodeMap & node_map{key_val.second};
		std::vector<NodeMap::Assignment> node_map_as_rel;
		node_map.as_relation(node_map_as_rel);
		NodeMap new_node_map(median.num_nodes, node_map.num_target_nodes());
		for (const auto & assignment : node_map_as_rel) {
			new_node_map.add_assignment(assignment.first, assignment.second);
		}
		new_node_map.add_assignment(median.num_nodes - 1, best_config.at(graph_id));
		node_map = new_node_map;
	}

	// Increase overall number of increases.
	num_increase_order_++;
}


}

#endif /* MEDIAN_SRC_MEDIAN_GRAPH_ESTIMATOR_IPP_ */
