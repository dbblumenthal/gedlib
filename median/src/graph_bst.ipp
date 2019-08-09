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
 * @file graph_bst.ipp
 * @brief ged::GraphBST class definition.
 */

#ifndef MEDIAN_SRC_GRAPH_BST_IPP_
#define MEDIAN_SRC_GRAPH_BST_IPP_


namespace ged {

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
GraphBST(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge) :
ged_env_{ged_env},
mge_{mge},
lower_bound_method_{Options::GEDMethod::BRANCH_FAST},
lower_bound_options_(""),
upper_bound_method_{Options::GEDMethod::IPFP},
upper_bound_options_(""),
max_cluster_size_{10},
focal_graphs_("MEDIANS"),
cutoff_{0.5},
print_to_stdout_{2},
root_{nullptr},
focal_graph_ids_(),
verified_graph_ids_(),
filtered_graph_ids_(),
undecided_graph_ids_(),
init_time_(),
query_time_(),
num_ged_evals_focal_graphs_{0},
num_ged_evals_data_graphs_{0} {
	if (ged_env_ == nullptr) {
		throw Error("The environment pointer passed to the constructor of ged::GraphBST is null.");
	}
	else if (not ged_env_->initialized()) {
		throw Error("The environment is uninitialized. Call ged::GEDEnv::init() before passing it to the constructor of ged::GraphBST.");
	}
	if (mge_ != nullptr) {
		if (ged_env_ != mge_->get_ged_env()) {
			throw Error("The environment pointers passed to the constructors of ged::GraphBST and ged::MedianGraphEstimator don't coincide.");
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
~GraphBST() {
	if (root_ != nullptr) {
		delete root_;
		root_ = nullptr;
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_default_options_() {
	focal_graphs_ = "MEDIANS";
	max_cluster_size_ = 10;
	cutoff_ = 0.5;
	print_to_stdout_ = 2;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_options(const std::string & options) {
	set_default_options_();
	std::map<std::string, std::string> options_map;
	util::options_string_to_options_map(options, options_map);
	for (const auto & option : options_map) {
		if (option.first == "focal-graphs") {
			focal_graphs_ = option.second;
			if (option.second != "MEDOIDS" and option.second != "MEDIANS" and option.second != "CENTERS") {
				throw ged::Error(std::string("Invalid argument ") + option.second + " for option focal-graphs. Usage: options = \"[--focal-graphs MEDIANS|MEDOIDS|CENTERS] [...]\"");
			}
		}
		else if (option.first == "max-cluster-size") {
			try {
				max_cluster_size_ = std::stoul(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option max-cluster-size. Usage: options = \"[--max-cluster-size <convertible to int greater 0>]\"");
			}
			if (max_cluster_size_ == 0) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option max-cluster-size. Usage: options = \"[--max-cluster-size <convertible to int greater 0>]\"");
			}
		}
		else if (option.first == "cutoff") {
			try {
				cutoff_ = std::stod(option.second);
			}
			catch (...) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option cutoff. Usage: options = \"[--cutoff <convertible to double greater 0 and smaller equal 1>]\"");
			}
			if (cutoff_ <= 0 or cutoff_ > 1) {
				throw Error(std::string("Invalid argument \"") + option.second + "\" for option cutoff. Usage: options = \"[--cutoff <convertible to double greater 0 and smaller equal 1>]\"");
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
			std::string valid_options("[--focal-graphs <arg>] [--max-cluster-size <arg>] [--cutoff <arg>] [--stdout <arg>]");
			throw Error(std::string("Invalid option \"") + option.first + "\". Usage: options = \"" + valid_options + "\"");
		}
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
init(std::vector<GEDGraph::GraphID> graph_ids, std::vector<GEDGraph::GraphID> focal_graph_ids) {

	// Start recording the initialization time.
	auto start = std::chrono::high_resolution_clock::now();

	// Unitialize the BST if it has already been initialized.
	if (root_ != nullptr) {
		delete root_;
		root_ = nullptr;
	}

	// Set the focal graph IDs and ensure that enough space has been allocated.
	focal_graph_ids_ = focal_graph_ids;
	if (focal_graph_ids_.size() < 2 * graph_ids.size()) {
		throw Error("Not enough focal graph IDs provide; ensure focal_graph_ids.size() >= 2 * graph_ids.size().");
	}

	// Initialize the root of the BST.
	if (print_to_stdout_ == 2) {
		std::cout << "Initializing the root: ... " << std::flush;
	}
	std::size_t pos_next_focal_graph_id{0};
	GEDGraph::GraphID root_id{focal_graph_ids_.at(pos_next_focal_graph_id++)};
	std::map<GEDGraph::GraphID, double> distances_from_root;
	compute_new_focal_graph_(graph_ids, root_id, distances_from_root);
	double radius{0};
	for (const auto & key_val : distances_from_root) {
		if (key_val.second > radius) {
			radius = key_val.second;
		}
	}
	root_ = new BSTNode_();
	init_bst_node_(root_, root_id, radius, graph_ids, distances_from_root);
	if (print_to_stdout_ == 2) {
		std::cout << "done.\n";
	}

	// Recursively build the BST and mark the BST as initialized.
	if (root_->cluster.size() > max_cluster_size_) {
		ProgressBar progress(graph_ids.size());
		if (print_to_stdout_ == 2) {
			std::cout << "\rInitializing BST: " << progress << std::flush;
		}
		split_(root_, pos_next_focal_graph_id, progress);
		std::cout << "\n";
	}

	// Set the initialization time.
	auto end = std::chrono::high_resolution_clock::now();
	init_time_ = start - end;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
init_bst_node_(BSTNode_ * bst_node, GEDGraph::GraphID focal_graph_id, double radius, const std::vector<GEDGraph::GraphID> & cluster, const std::map<GEDGraph::GraphID, double> & distances_from_focal_graph) {
	bst_node->focal_graph_id = focal_graph_id;
	bst_node->radius = radius;
	bst_node->cluster = cluster;
	bst_node->distances_from_focal_graph = distances_from_focal_graph;
	bst_node->child_1 = nullptr;
	bst_node->child_2 = nullptr;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
compute_new_focal_graph_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID new_focal_graph_id, std::map<GEDGraph::GraphID, double> & distances_from_new_focal_graph) {
	distances_from_new_focal_graph.clear();
	if (focal_graphs_ == "MEDIANS") {
		// Determine the new focal graph.
		mge_->run(graph_ids, new_focal_graph_id);

		// Save the distances from the new focal graph.
		for (GEDGraph::GraphID graph_id : graph_ids) {
			distances_from_new_focal_graph[graph_id] = mge_->get_distance_from_median(graph_id);
		}
	}
	else {
		// Determine the new focal graph.
		ged_env_->set_method(upper_bound_method_, upper_bound_options_);
		GEDGraph::GraphID original_focal_graph_id{undefined()};
		double best_aggr_dist{std::numeric_limits<double>::infinity()};
		for (auto g_id : graph_ids) {
			double aggr_dist{0};
			for (auto h_id : graph_ids) {
				ged_env_->run_method(g_id, h_id);
				if (focal_graphs_ == "MEDOIDS")  {
					aggr_dist += ged_env_->get_upper_bound(g_id, h_id);
				}
				else {
					aggr_dist = std::max(aggr_dist, ged_env_->get_upper_bound(g_id, h_id));
				}
			}
			if (aggr_dist < best_aggr_dist) {
				best_aggr_dist = aggr_dist;
				original_focal_graph_id = g_id;
			}
		}

		// Save the distances from the new focal graph.
		for (GEDGraph::GraphID graph_id : graph_ids) {
			distances_from_new_focal_graph[graph_id] = ged_env_->get_upper_bound(original_focal_graph_id, graph_id);
		}

		// Load the new focal graph into the environment.
		ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> focal_graph(ged_env_->get_graph(original_focal_graph_id));
		ged_env_->load_exchange_graph(focal_graph, new_focal_graph_id);
		ged_env_->init(ged_env_->get_init_type());
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
split_(BSTNode_ * node, std::size_t pos_next_focal_graph_id, ProgressBar & progress) {

	// Get the IDs of the children's focal graphs.
	GEDGraph::GraphID child_id_2{focal_graph_ids_.at(pos_next_focal_graph_id++)};

	// Collect IDs of the graphs in the cluster that are far away from the current median.
	double min_dist{std::numeric_limits<double>::infinity()};
	for (const auto & key_val : node->distances_from_focal_graph) {
		min_dist = std::min(min_dist, key_val.second);
	}
	double cutoff{min_dist + cutoff_ * (node->radius - min_dist)};
	std::vector<GEDGraph::GraphID> cluster_2;
	std::map<GEDGraph::GraphID, bool> is_far_away;
	for (const auto & key_val : node->distances_from_focal_graph) {
		if (key_val.second >= cutoff) {
			cluster_2.emplace_back(key_val.first);
			is_far_away[key_val.first] = true;
		}
		else {
			is_far_away[key_val.first] = false;
		}
	}

	// Compute the new focal graph and its distances from the nodes in the cluster.
	std::map<GEDGraph::GraphID, double> distances_from_new_focal_graph;
	compute_new_focal_graph_(cluster_2, child_id_2, distances_from_new_focal_graph);
	for (GEDGraph::GraphID graph_id : node->cluster) {
		if (not is_far_away.at(graph_id)) {
			if (focal_graphs_ == "MEDIANS") {
				distances_from_new_focal_graph[graph_id] = mge_->compute_node_map_from_median(graph_id).induced_cost();
			}
			else {
				ged_env_->set_method(upper_bound_method_, upper_bound_options_);
				ged_env_->run_method(child_id_2, graph_id);
				distances_from_new_focal_graph[graph_id] = ged_env_->get_upper_bound(child_id_2, graph_id);
			}
		}
	}

	// Compute the new clusters and the new radii.
	std::vector<GEDGraph::GraphID> cluster_1;
	cluster_2.clear();
	std::map<GEDGraph::GraphID, double> distances_from_focal_graph_1;
	std::map<GEDGraph::GraphID, double> distances_from_focal_graph_2;
	double radius_1{0};
	double radius_2{0};
	for (GEDGraph::GraphID graph_id : node->cluster) {
		double distance_1{node->distances_from_focal_graph.at(graph_id)};
		double distance_2{distances_from_new_focal_graph.at(graph_id)};
		if (distance_1 <= distance_2) {
			cluster_1.emplace_back(graph_id);
			distances_from_focal_graph_1.emplace(graph_id, distance_1);
			radius_1 = std::max(radius_1, distance_1);
		}
		else {
			cluster_2.emplace_back(graph_id);
			distances_from_focal_graph_2.emplace(graph_id, distance_2);
			radius_2 = std::max(radius_2, distance_2);
		}
	}


	// Initialize the children.
	node->child_1 = new BSTNode_();
	init_bst_node_(node->child_1, node->focal_graph_id, radius_1, cluster_1, distances_from_focal_graph_1);
	node->child_2 = new BSTNode_();
	init_bst_node_(node->child_2, child_id_2, radius_2, cluster_2, distances_from_focal_graph_2);

	// Print information.
	if (print_to_stdout_ == 2) {
		progress.increment();
		std::cout << "\rInitializing BST: " << progress << std::flush;
	}

	// Recursively call split_() on the children.
	if (node->child_1->cluster.size() > max_cluster_size_) {
		pos_next_focal_graph_id = split_(node->child_1, pos_next_focal_graph_id, progress);
	}
	if (node->child_2->cluster.size() > max_cluster_size_) {
		pos_next_focal_graph_id = split_(node->child_2, pos_next_focal_graph_id, progress);
	}

	// Return the position of the next focal graph ID.
	return pos_next_focal_graph_id;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
process_range_query(GEDGraph::GraphID query_graph_id, double threshold) {

	if (print_to_stdout_ > 0) {
		std::cout << "\n===========================================================\n";
		std::cout << "Processing range query.\n";
		std::cout << "-----------------------------------------------------------\n";
	}

	// Ensure that the BST has been initialized.
	if (root_ == nullptr) {
		throw Error("The bisector tree has not been initialized. Call init() or load() before calling process_range_query().");
	}

	// Start recording the query time.
	auto start = std::chrono::high_resolution_clock::now();

	// Clear the result variables.
	verified_graph_ids_.clear();
	filtered_graph_ids_.clear();
	undecided_graph_ids_.clear();
	num_ged_evals_focal_graphs_ = 0;
	num_ged_evals_data_graphs_ = 0;

	// Recursively process the query.
	check_(query_graph_id, threshold, root_, -1, -1);

	// Set the query time.
	auto end = std::chrono::high_resolution_clock::now();
	query_time_ = end - start;

	if (print_to_stdout_ > 0) {
		std::cout << "Threshold: " << threshold << "\n";
		std::cout << "Number of GED evaluations (total): " << num_ged_evals_focal_graphs_ + num_ged_evals_data_graphs_ << "\n";
		std::cout << "Number of GED evaluations (inner): " << num_ged_evals_focal_graphs_ << "\n";
		std::cout << "Number of GED evaluations (leafs): " << num_ged_evals_data_graphs_ << "\n";
		std::cout << "Number of verified data graphs: " << verified_graph_ids_.size() << "\n";
		std::cout << "Number of filtered data graphs: " << filtered_graph_ids_.size() << "\n";
		std::cout << "Number of undecided data graphs: " << undecided_graph_ids_.size() << "\n";
		std::cout << "Initialization time: " << get_init_time() << "\n";
		std::cout << "Query time: " << get_query_time() << "\n";
		std::cout << "===========================================================\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
check_(GEDGraph::GraphID query_graph_id, double threshold, const BSTNode_ * node, double lower_bound, double upper_bound) {

	// Update number of GED evaluations.
	bool recompute_bounds{false};
	if (lower_bound < 0) {
		num_ged_evals_focal_graphs_++;
		recompute_bounds = true;
	}

	// Compute lower bound.
	if (recompute_bounds) {
		ged_env_->set_method(lower_bound_method_, lower_bound_options_);
		ged_env_->run_method(query_graph_id, node->focal_graph_id);
		lower_bound = ged_env_->get_lower_bound(query_graph_id, node->focal_graph_id);
	}

	// Use lower bound for filtering.
	if (lower_bound > threshold + node->radius) {
		for (GEDGraph::GraphID graph_id : node->cluster) {
			filtered_graph_ids_.emplace_back(graph_id);
		}
		return;
	}

	// Compute upper bound if the upper bound method is different from the lower bound method.
	if (recompute_bounds) {
		if (lower_bound_method_ != upper_bound_method_ or lower_bound_options_ != upper_bound_options_) {
			ged_env_->set_method(upper_bound_method_, upper_bound_options_);
			ged_env_->run_method(query_graph_id, node->focal_graph_id);
		}
		upper_bound = ged_env_->get_upper_bound(query_graph_id, node->focal_graph_id);
	}

	// Use upper bound for verifying.
	if (upper_bound <= threshold - node->radius) {
		for (GEDGraph::GraphID graph_id : node->cluster) {
			verified_graph_ids_.emplace_back(graph_id);
		}
		return;
	}

	// The current node is a leaf that has neither been filtered nor verified.
	if (node->child_1 == nullptr) {
		for (GEDGraph::GraphID graph_id : node->cluster) {
			if (lower_bound > threshold + node->distances_from_focal_graph.at(graph_id)) {
				filtered_graph_ids_.emplace_back(graph_id);
				continue;
			}
			if (upper_bound <= threshold - node->distances_from_focal_graph.at(graph_id)) {
				verified_graph_ids_.emplace_back(graph_id);
				continue;
			}
			num_ged_evals_data_graphs_++;
			ged_env_->set_method(lower_bound_method_, lower_bound_options_);
			ged_env_->run_method(query_graph_id, graph_id);
			if (ged_env_->get_lower_bound(query_graph_id, graph_id) > threshold) {
				filtered_graph_ids_.emplace_back(graph_id);
			}
			else {
				if (lower_bound_method_ != upper_bound_method_ or lower_bound_options_ != upper_bound_options_) {
					ged_env_->set_method(upper_bound_method_, upper_bound_options_);
					ged_env_->run_method(query_graph_id, graph_id);
				}
				if (ged_env_->get_upper_bound(query_graph_id, graph_id) <= threshold) {
					verified_graph_ids_.emplace_back(graph_id);
				}
				else {
					undecided_graph_ids_.emplace_back(graph_id);
				}
			}
		}
		return;
	}

	// The current node is an inner node that has neither been filtered nor verified.
	check_(query_graph_id, threshold, node->child_1, lower_bound, upper_bound);
	check_(query_graph_id, threshold, node->child_2, -1, -1);

}

template<>
void
GraphBST<GXLNodeID, GXLLabel, GXLLabel>::
serialize_(const BSTNode_ * node, const std::string & node_id, std::size_t max_cluster_size, const std::string & focal_graph_dir, std::map<std::string, std::string> & config) const {

	// Save the focal graph of the current node.
	if (focal_graph_dir != "") {
		ged_env_->save_as_gxl_graph(node->focal_graph_id, focal_graph_dir + "/" + ged_env_->get_graph_name(node->focal_graph_id));
	}

	// Save the radius and the cluster in the configuration map.
	config[node_id] = ged_env_->get_graph_name(node->focal_graph_id) + ";" + std::to_string(node->radius);
	for (GEDGraph::GraphID graph_id : node->cluster) {
		config[node_id] += ";" + ged_env_->get_graph_name(graph_id);
		double distance_from_focal_graph{node->distances_from_focal_graph.at(graph_id)};
		config[node_id] += "," + std::to_string(distance_from_focal_graph);
	}

	if (node->cluster.size() <= max_cluster_size) {
		config[node_id + "1"] = "NULL";
		config[node_id + "2"] = "NULL";
	}
	else {
		serialize_(node->child_1, node_id + "1", max_cluster_size, focal_graph_dir, config);
		serialize_(node->child_2, node_id + "2", max_cluster_size, focal_graph_dir, config);
	}
}

template<>
void
GraphBST<GXLNodeID, GXLLabel, GXLLabel>::
save(const std::string & bst_file_name, const std::string & focal_graph_dir, std::size_t max_cluster_size) const {
	// Sanity check.
	if (root_ == nullptr) {
		throw Error("The bisector tree has not been initialized. Call init() or load() before calling save().");
	}

	// Serialize the tree.
	std::map<std::string, std::string> config;
	serialize_(root_, "0", std::max(max_cluster_size_, max_cluster_size), focal_graph_dir, config);

	// Write the configuration file.
	util::save_as_config_file(bst_file_name, config);
}

template<>
void
GraphBST<GXLNodeID, GXLLabel, GXLLabel>::
de_serialize_(const std::string & node_id, const std::map<std::string, std::string> & config, const std::string & data_graph_dir, const std::string & focal_graph_dir,
		Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
		const std::unordered_set<std::string> & irrelevant_node_attributes, const std::unordered_set<std::string> & irrelevant_edge_attributes,
		BSTNode_ * & node, std::map<std::string, GEDGraph::GraphID> & graph_name_to_id) {

	// Return if the entry in the configuration file is NULL.
	if (config.at(node_id) == "NULL") {
		return;
	}

	// Tokenize the information from the configuration file.
	std::vector<std::string> node_info;
	util::tokenize(config.at(node_id), ';', node_info);

	// Load the cluster and load the data graphs into the environment, if the current node is the root.
	std::vector<GEDGraph::GraphID> cluster;
	std::map<GEDGraph::GraphID, double> distances_from_focal_graph;
	for (std::size_t pos{2}; pos < node_info.size(); pos++) {
		std::vector<std::string> gxl_dist;
		util::tokenize(node_info.at(pos), ',', gxl_dist);
		GEDGraph::GraphID graph_id{undefined()};
		if (node_id == "0") {
			graph_id = ged_env_->load_gxl_graph(data_graph_dir, gxl_dist.at(0), node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, undefined(), "no_class");
			graph_name_to_id[gxl_dist.at(0)] = graph_id;
		}
		else {
			graph_id = graph_name_to_id.at(gxl_dist.at(0));
		}
		cluster.emplace_back(graph_id);
		distances_from_focal_graph.emplace(graph_id, std::stod(gxl_dist.at(1)));
	}

	// Load the focal graph into the environment.
	GEDGraph::GraphID focal_graph_id{ged_env_->load_gxl_graph(focal_graph_dir, node_info.at(0), node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, undefined(), "no_class")};

	// Get the radius.
	double radius{std::stod(node_info.at(1))};

	// Initialize the current node.
	node = new BSTNode_();
	init_bst_node_(node, focal_graph_id, radius, cluster, distances_from_focal_graph);

	// Continue with the children.
	de_serialize_(node_id + "1", config, data_graph_dir, focal_graph_dir, node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, node->child_1, graph_name_to_id);
	de_serialize_(node_id + "2", config, data_graph_dir, focal_graph_dir, node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, node->child_2, graph_name_to_id);
}

template<>
void
GraphBST<GXLNodeID, GXLLabel, GXLLabel>::
load(const std::string & bst_file_name, const std::string & data_graph_dir, const std::string & focal_graph_dir,
		Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
		const std::unordered_set<std::string> & irrelevant_node_attributes, const std::unordered_set<std::string> & irrelevant_edge_attributes) {

	if (print_to_stdout_ > 0) {
		std::cout << "\n===========================================================\n";
		std::cout << "Loading bisector tree from file.\n";
		std::cout << "-----------------------------------------------------------\n";
	}

	// Start recording the initialization time.
	auto start = std::chrono::high_resolution_clock::now();

	// Unitialize the BST if it has already been initialized.
	if (root_ != nullptr) {
		delete root_;
		root_ = nullptr;
	}

	// Parse the configuration file.
	std::map<std::string, std::string> config;
	util::parse_config_file(bst_file_name, config);

	// De-serialize the tree.
	std::map<std::string, GEDGraph::GraphID> graph_name_to_id;
	if (print_to_stdout_ == 2) {
		std::cout << "De-serializing the tree: ... " << std::flush;
	}
	de_serialize_("0", config, data_graph_dir, focal_graph_dir, node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, root_, graph_name_to_id);
	if (print_to_stdout_ == 2) {
		std::cout << "done.\n";
	}
	ged_env_->init(ged_env_->get_init_type(), (print_to_stdout_ == 2));

	// Set the initialization time.
	auto end = std::chrono::high_resolution_clock::now();
	init_time_ = end - start;
	if (print_to_stdout_ > 0) {
		std::cout << "===========================================================\n";
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_lower_bound_method(Options::GEDMethod lower_bound_method, const std::string & lower_bound_options) {
	lower_bound_method_ = lower_bound_method;
	lower_bound_options_ = lower_bound_options;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
set_upper_bound_method(Options::GEDMethod upper_bound_method, const std::string & upper_bound_options) {
	upper_bound_method_ = upper_bound_method;
	upper_bound_options_ = upper_bound_options;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const std::vector<GEDGraph::GraphID> &
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_verified_graph_ids() const {
	return verified_graph_ids_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const std::vector<GEDGraph::GraphID> &
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_filtered_graph_ids() const {
	return filtered_graph_ids_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
const std::vector<GEDGraph::GraphID> &
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_undecided_graph_ids() const {
	return undecided_graph_ids_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_init_time() const {
	return init_time_.count();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
double
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_query_time() const {
	return query_time_.count();
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_ged_evals() const {
	return num_ged_evals_focal_graphs_ + num_ged_evals_data_graphs_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_ged_evals_focal_graphs() const {
	return num_ged_evals_focal_graphs_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
get_num_ged_evals_data_graphs() const {
	return num_ged_evals_data_graphs_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
BSTNode_::
BSTNode_() :
focal_graph_id{undefined()},
radius{std::numeric_limits<double>::infinity()},
cluster(),
distances_from_focal_graph(),
child_1{nullptr},
child_2{nullptr} {}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
BSTNode_::
~BSTNode_() {
	if (child_1 != nullptr) {
		delete child_1;
		child_1 = nullptr;
	}
	if (child_2 != nullptr) {
		delete child_2;
		child_2 = nullptr;
	}
}


}


#endif /* MEDIAN_SRC_GRAPH_BST_IPP_ */
