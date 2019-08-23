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
 * @file graph_clustering_heuristic.ipp
 * @brief ged::GraphClusteringHeuristic class definition.
 */

#ifndef MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTIC_IPP_
#define MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTIC_IPP_

namespace ged {

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
GraphClusteringHeuristic(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge):
ged_env_{ged_env},
mge_{mge},
ged_method_{Options::GEDMethod::BRANCH_UNIFORM},
ged_options_(""),
clustering_method_("K-MEDIANS"),
init_type_("K-MEANS++"),
num_random_inits_{10},
use_real_randomness_{true},
seed_{0},
time_limit_in_sec_{0},
epsilon_{0.0001},
max_itrs_{100},
print_to_stdout_{2},
assigned_focal_graph_ids_(),
node_maps_from_assigned_focal_graphs_(),
cluster_sums_of_distances_(),
cluster_radii_(),
sum_of_distances_{0},
runtime_(),
itrs_() {
	if (ged_env_ == nullptr) {
		throw Error("The environment pointer passed to the constructor of ged::GraphClusteringHeuristic is null.");
	}
	else if (not ged_env_->initialized()) {
		throw Error("The environment is uninitialized. Call ged::GEDEnv::init() before passing it to the constructor of ged::GraphClusteringHeuristic.");
	}
	else if (mge_ != nullptr) {
		if (ged_env_ != mge_->get_ged_env()) {
			throw Error("The pointers passed to the constructors of ged::GraphClusteringHeuristic and ged::MedianGraphEstimators don't coincide.");
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_options(const std::string & options) {
	set_default_options_();
	std::map<std::string, std::string> options_map;
	util::options_string_to_options_map(options, options_map);
	for (const auto & option : options_map) {
		if (option.first == "clustering-method") {
			clustering_method_ = option.second;
			if (option.second != "K-MEDIANS" and option.second != "K-MEDOIDS") {
				throw ged::Error(std::string("Invalid argument ") + option.second + " for option clustering-method. Usage: options = \"[--clustering-method K-MEDIANS|K-MEDOIDS] [...]\"");
			}
		}
		else if (option.first == "init-type") {
			init_type_ = option.second;
			if (option.second != "K-MEANS++" and option.second != "CLUSTERS") {
				throw ged::Error(std::string("Invalid argument ") + option.second + " for option init-type. Usage: options = \"[--init-type K-MEANS++|CLUSTERS] [...]\"");
			}
		}
		else if (option.first == "random-inits") {
			try {
				num_random_inits_ = std::stoul(option.second);
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
		else if (option.first == "seed") {
			try {
				seed_ = std::stoul(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option seed. Usage: options = \"[--seed <convertible to int greater equal 0>] [...]");
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
		else if (option.first == "max-itrs") {
			try {
				max_itrs_ = std::stoi(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option max-itrs. Usage: options = \"[--max-itrs <convertible to int>] [...]");
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
		else {
			std::string valid_options("[--clustering-method <arg>] [--init-type <arg>] ");
			valid_options += "[--random-inits <arg>] [--randomness <arg>] [--seed <arg>] [--stdout <arg>] ";
			valid_options += "[--time-limit <arg>] [--max-itrs <arg>] [--epsilon <arg>] ";
			throw Error(std::string("Invalid option \"") + option.first + "\". Usage: options = \"" + valid_options + "\"");
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_ged_method(Options::GEDMethod ged_method, const std::string ged_options) {
	ged_method_ = ged_method;
	ged_options_ = ged_options;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
run(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids) {

	// Sanity checks.
	if (graph_ids.empty()) {
		throw Error("Empty vector of graph IDs, unable to compute clusters.");
	}
	if (focal_graph_ids.size() < 2) {
		throw Error("Vector of focal graph IDs has size < 2, unable to compute clusters.");
	}
	if ((mge_ == nullptr) and (clustering_method_ == "K-MEDIANS") and (max_itrs_ > 0)) {
		throw Error("No median graph estimator availabel, unable to compute clusters.");
	}

	// Start timer and record start time.
	auto start = std::chrono::high_resolution_clock::now();
	Timer timer(time_limit_in_sec_);

	// Reset information about iterations and number of times the median decreases and increases.
	itrs_ = std::vector<std::size_t>(num_random_inits_, 0);

	// Seed the random number generator.
	std::mt19937 urng;
	if (use_real_randomness_) {
		std::random_device rng;
		urng.seed(rng());
	}
	else {
		urng.seed(seed_);
	}

	// Select the employed GED method.
	ged_env_->set_method(ged_method_, ged_options_);

	// Initialize the best clustering.
	double best_sum_of_distances{std::numeric_limits<double>::infinity()};
	std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> best_focal_graphs;
	std::map<GEDGraph::GraphID, GEDGraph::GraphID> best_assigned_focal_graph_ids;
	std::map<GEDGraph::GraphID, NodeMap> best_node_maps_from_assigned_focal_graphs;
	std::map<GEDGraph::GraphID, double> best_distances_from_assigned_focal_graphs;
	std::map<GEDGraph::GraphID, double> best_cluster_sums_of_distances;

	// Cluster the graphs.
	for (std::size_t random_init{0}; random_init < num_random_inits_; random_init++) {

		// Generate the initial focal graphs and the initial clusters.
		std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> focal_graphs;
		compute_initial_clusters_(graph_ids, focal_graph_ids, urng, focal_graphs);

		// Run Lloyd's algorithm.
		bool converged{false};
		for (; not termination_criterion_met_(converged, timer, itrs_.at(random_init)); itrs_[random_init]) {
			converged = update_clusters_(graph_ids, focal_graph_ids);
			if (not converged) {
				update_focal_graphs_(focal_graph_ids, focal_graphs);
			}
		}

		// Update the best clustering.
		if (sum_of_distances_ < best_sum_of_distances - epsilon_) {
			best_focal_graphs = focal_graphs;
			best_sum_of_distances = sum_of_distances_;
			best_assigned_focal_graph_ids = assigned_focal_graph_ids_;
			best_node_maps_from_assigned_focal_graphs = node_maps_from_assigned_focal_graphs_;
			best_cluster_sums_of_distances = cluster_sums_of_distances_;
		}
	}

	// Store the best encountered clustering.
	for (const auto & key_val : best_focal_graphs) {
		ged_env_->load_exchange_graph(key_val.second, key_val.first);
	}
	ged_env_->init(ged_env_->get_init_type());
	sum_of_distances_ = best_sum_of_distances;
	assigned_focal_graph_ids_ = best_assigned_focal_graph_ids;
	node_maps_from_assigned_focal_graphs_ = best_node_maps_from_assigned_focal_graphs;
	cluster_sums_of_distances_ = best_cluster_sums_of_distances;

	// Compute the cluster radii.
	compute_cluster_radii_(focal_graph_ids);

	// Record the runtime.
	auto end = std::chrono::high_resolution_clock::now();
	runtime_ = start - end;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_default_options_() {
	clustering_method_ = "K-MEDIANS";
	init_type_ = "K-MEANS++";
	num_random_inits_ = 10;
	use_real_randomness_ = true;
	seed_ = 0;
	time_limit_in_sec_ = 0;
	epsilon_ = 0.0001;
	max_itrs_ = 100;
	print_to_stdout_ = 2;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_initial_clusters_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
		std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs) {

	// Generate the focal graphs.
	focal_graphs.clear();
	if (init_type_ == "K-MEANS++") {
		compute_initial_focal_graphs_k_means_plus_plus_(graph_ids, focal_graph_ids, urng, focal_graphs);
	}
	else {
		compute_initial_focal_graphs_cluster_sampling_(graph_ids, focal_graph_ids, urng, focal_graphs);
	}

	// Initialize the clusters.
	assigned_focal_graph_ids_.clear();
	node_maps_from_assigned_focal_graphs_.clear();
	cluster_sums_of_distances_.clear();
	NodeMap empty_node_map(0, 0);
	for (GEDGraph::GraphID graph_id : graph_ids) {
		assigned_focal_graph_ids_.emplace(graph_ids, undefined());
		node_maps_from_assigned_focal_graphs_.emplace(graph_id, empty_node_map);
	}
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
		cluster_sums_of_distances_.emplace(focal_graph_id, std::numeric_limits<double>::infinity());
		cluster_radii_.emplace(focal_graph_id, std::numeric_limits<double>::infinity());
	}

}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_initial_focal_graphs_k_means_plus_plus_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
		std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs) {

	// IDs of the graphs that will be selected as focal graphs.
	std::vector<GEDGraph::GraphID> original_focal_graph_ids;

	// Flags that indicate for all graphs whether they have been selected as focal graphs.
	std::vector<bool> already_selected(graph_ids.size(), false);

	// Select the first focal graph.
	std::uniform_int_distribution<std::size_t> unif_dist(0, graph_ids.size() - 1);
	GEDGraph::GraphID selected_graph_id_pos{unif_dist(urng)};
	original_focal_graph_ids.emplace_back(graph_ids.at(selected_graph_id_pos));
	focal_graphs.emplace(focal_graph_ids.at(focal_graphs.size()), ged_env_->get_graph(original_focal_graph_ids.back()));
	already_selected[selected_graph_id_pos] = true;

	// Select all remaining focal graphs.
	while (focal_graphs.size()  < focal_graph_ids.size()) {

		// Compute the weights that yield the probability distribution.
		std::vector<double> weights(graph_ids.size(), std::numeric_limits<double>::infinity());
		for (std::size_t graph_id_pos{0}; graph_id_pos < graph_ids.size(); graph_id_pos++) {
			if (already_selected.at(graph_id_pos)) {
				weights[graph_id_pos] = 0;
				continue;
			}
			for (const auto & focal_graph_id : original_focal_graph_ids) {
				ged_env_->run_method(focal_graph_id, graph_ids.at(graph_id_pos));
				weights[graph_id_pos] = std::min(weights.at(graph_id_pos), ged_env_->get_upper_bound(focal_graph_id, graph_ids.at(graph_id_pos)));
			}
		}

		// Draw the next focal graph from the obtained distribution.
		std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
		selected_graph_id_pos = dist(urng);
		original_focal_graph_ids.emplace_back(graph_ids.at(selected_graph_id_pos));
		focal_graphs.emplace(focal_graph_ids.at(focal_graphs.size()), ged_env_->get_graph(original_focal_graph_ids.back()));
		already_selected[selected_graph_id_pos] = true;
	}

	// Load the focal graphs into the environment.
	for (std::size_t cluster_id{0}; cluster_id < focal_graph_ids.size(); cluster_id++) {
		ged_env_->load_exchange_graph(focal_graphs.at(cluster_id), focal_graph_ids.at(cluster_id));
	}
	ged_env_->init(ged_env_->get_init_type());
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_initial_focal_graphs_cluster_sampling_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
		std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs) {

	// Determine the cluster size.
	std::size_t cluster_size{graph_ids.size() / focal_graph_ids.size()};

	// Shuffle the graph IDs.
	std::vector<GEDGraph::GraphID> shuffled_graph_ids(graph_ids);
	std::shuffle(shuffled_graph_ids.begin(), shuffled_graph_ids.end(), urng);

	// Initialize the clustering.
	std::vector<std::vector<GEDGraph::GraphID>> clustering;

	// Populate the clustering.
	std::size_t pos{0};
	while (clustering.size()  < focal_graph_ids.size() - 1) {
		clustering.emplace_back();
		while (pos++ < (clustering.size() + 1) * cluster_size) {
			clustering.back().emplace_back(shuffled_graph_ids.at(pos));
		}
	}
	clustering.emplace_back();
	while (pos++ < shuffled_graph_ids.size()) {
		clustering.back().emplace_back(shuffled_graph_ids.at(pos));
	}

	// Compute the initial focal graphs as medians or medoids of the initial clusters.
	for (std::size_t cluster_id{0}; cluster_id < clustering.size(); cluster_id++) {
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> focal_graph;
		std::map<GEDGraph::GraphID, NodeMap> new_node_maps_from_focal_graph;
		if (clustering_method_ == "K-MEDIANS") {
			compute_median_(clustering.at(cluster_id), focal_graph_ids.at(cluster_id), focal_graph, new_node_maps_from_focal_graph);
		}
		else {
			compute_medoid_(clustering.at(cluster_id), focal_graph_ids.at(cluster_id), focal_graph, new_node_maps_from_focal_graph);
		}
		focal_graphs.emplace(focal_graph_ids.at(cluster_id), focal_graph);
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
termination_criterion_met_(bool converged, const Timer & timer, std::size_t itr) const {
	return converged or timer.expired() or (max_itrs_ >= 0 ? itr >= max_itrs_ : false);
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
bool
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_clusters_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids) {

	// Re-select the employed GED method if K-MEDIANS is used (the median graph estimator might have used another method).
	if (clustering_method_ == "K-MEDIANS") {
		ged_env_->set_method(ged_method_, ged_options_);
	}

	// Determine the closest focal graph for each graph.
	bool clusters_modified{false};
	for (GEDGraph::GraphID graph_id : graph_ids) {
		GEDGraph::GraphID old_assigned_focal_graph_id{assigned_focal_graph_ids_.at(graph_id)};
		for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
			ged_env_->run_method(focal_graph_id, graph_id);
			if (ged_env_->get_upper_bound(focal_graph_id, graph_id) < node_maps_from_assigned_focal_graphs_.at(graph_id).induced_cost() - epsilon_) {
				assigned_focal_graph_ids_[graph_id] = focal_graph_id;
				node_maps_from_assigned_focal_graphs_[graph_id] = ged_env_->get_node_map(focal_graph_id, graph_id);
			}
		}
		if (assigned_focal_graph_ids_.at(graph_id) != old_assigned_focal_graph_id) {
			clusters_modified = true;
		}
	}

	// Update the sums of distances for the clusters and the overall sum of distances.
	if (clusters_modified) {
		sum_of_distances_ = 0;
		for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
			cluster_sums_of_distances_.at(focal_graph_id) = 0;
		}
		for (const auto & key_val : assigned_focal_graph_ids_) {
			double new_cost{node_maps_from_assigned_focal_graphs_.at(key_val.first).induced_cost()};
			cluster_sums_of_distances_.at(key_val.second) += new_cost;
			sum_of_distances_ += new_cost;
		}
	}


	// Return true if the clusters were modified.
	return clusters_modified;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
update_focal_graphs_(const std::vector<GEDGraph::GraphID> & focal_graph_ids, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs) {

	// Collect the clusters.
	std::map<GEDGraph::GraphID, std::vector<GEDGraph::GraphID>> clustering;
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
		clustering.emplace(focal_graph_id, std::vector<GEDGraph::GraphID>());
	}
	for (const auto & key_val : assigned_focal_graph_ids_) {
		clustering.at(key_val.second).emplace_back(key_val.first);
	}

	// Update the focal graphs, the sums of distances for the clusters, and the node maps.
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
		std::map<GEDGraph::GraphID, NodeMap> new_node_maps_from_focal_graph;
		if (clustering_method_ == "K-MEDIANS") {
			cluster_sums_of_distances_.at(focal_graph_id) = compute_median_(clustering.at(focal_graph_id), focal_graph_id, focal_graphs.at(focal_graph_id), new_node_maps_from_focal_graph);
		}
		else {
			cluster_sums_of_distances_.at(focal_graph_id) = compute_medoid_(clustering.at(focal_graph_id), focal_graph_id, focal_graphs.at(focal_graph_id), new_node_maps_from_focal_graph);
		}
		for (GEDGraph::GraphID graph_id : clustering.at(focal_graph_id)) {
			node_maps_from_assigned_focal_graphs_.at(graph_id) = new_node_maps_from_focal_graph.at(graph_id);
		}
	}

	// Update the overall sum of distances.
	sum_of_distances_ = 0;
	for (const auto & key_val : cluster_sums_of_distances_) {
		sum_of_distances_ += key_val.second;
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_median_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID median_id,
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median, std::map<GEDGraph::GraphID, NodeMap> & node_maps_from_median) const {

	// Determine the median and load it into the environment.
	mge_->run(graph_ids, median_id);

	// Update the median's ExchangeGraph representation.
	median = ged_env_->get_graph(median_id);

	// Update the node maps within the cluster.
	node_maps_from_median.clear();
	for (GEDGraph::GraphID graph_id : graph_ids) {
		node_maps_from_median.emplace(graph_id, mge_->get_node_map_from_median(graph_id));
	}

	// Return the sum of distances.
	return mge_->get_sum_of_distances();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_medoid_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID medoid_id,
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & medoid, std::map<GEDGraph::GraphID, NodeMap> & node_maps_from_medoid) const {

	// Determine the medoid.
	GEDGraph::GraphID original_medoid_id{undefined()};
	double best_sum_of_distances{std::numeric_limits<double>::infinity()};
	for (auto g_id : graph_ids) {
		double sum_of_distances{0};
		for (auto h_id : graph_ids) {
			ged_env_->run_method(g_id, h_id);
			sum_of_distances += ged_env_->get_upper_bound(g_id, h_id);
		}
		if (sum_of_distances < best_sum_of_distances) {
			best_sum_of_distances = sum_of_distances;
			original_medoid_id = g_id;
		}
	}

	// Update the medoid's exchange graph representation.
	medoid = ged_env_->get_graph(original_medoid_id);

	// Update the node maps within the cluster.
	node_maps_from_medoid.clear();
	for (auto graph_id : graph_ids) {
		node_maps_from_medoid.emplace(graph_id, ged_env_->get_node_map(original_medoid_id, graph_id));
	}

	// Load the medoid into the environment.
	ged_env_->load_exchange_graph(medoid, medoid_id);
	ged_env_->init(ged_env_->get_init_type());

	// Return the sum of distances.
	return best_sum_of_distances;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_cluster_radii_(const std::vector<GEDGraph::GraphID> & focal_graph_ids) {

	// Initialize the cluster radii.
	for (GEDGraph::GraphID focal_graph_id : focal_graph_ids) {
		cluster_radii_.emplace(focal_graph_id, 0);
	}

	// Compute the radii.
	for (const auto & key_val : assigned_focal_graph_ids_) {
		double old_radius{cluster_radii_.at(key_val.second)};
		cluster_radii_.at(key_val.second) = std::max(old_radius, node_maps_from_assigned_focal_graphs_.at(key_val.first).induced_cost());
	}
}

}


#endif /* MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTIC_IPP_ */
