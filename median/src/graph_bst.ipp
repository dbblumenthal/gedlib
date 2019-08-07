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
GraphBST(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel> * clustering_heuristic) :
ged_env_{ged_env},
clustering_heuristic_{clustering_heuristic},
lower_bound_method_{Options::GEDMethod::BRANCH_FAST},
lower_bound_options_(""),
upper_bound_method_{Options::GEDMethod::IPFP},
upper_bound_options_(""),
root_{nullptr},
focal_graph_ids_(),
verified_graph_ids_(),
filtered_graph_ids_(),
undecided_graph_ids_(),
init_time_(),
query_time_(),
num_ged_evals_{0} {
	if (ged_env_ == nullptr) {
		throw Error("The environment pointer passed to the constructor of ged::GraphBST is null.");
	}
	else if (not ged_env_->initialized()) {
		throw Error("The environment is uninitialized. Call ged::GEDEnv::init() before passing it to the constructor of ged::GraphBST.");
	}
	if (clustering_heuristic_ == nullptr) {
		throw Error("The clustering heuristic pointer passed to the constructor of ged::GraphBST is null.");
	}
	if (ged_env_ != clustering_heuristic_->get_ged_env()) {
		throw Error("The environment pointers passed to the constructors of ged::GraphBST and ged::GraphClusteringHeuristic don't coincide.");
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
init(std::vector<GEDGraph::GraphID> graph_ids, std::vector<GEDGraph::GraphID> focal_graph_ids, std::size_t max_cluster_size) {

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
	std::size_t pos_next_focal_graph_id{0};
	GEDGraph::GraphID root_id{focal_graph_ids_.at(pos_next_focal_graph_id++)};
	clustering_heuristic_->run(graph_ids, {root_id});
	root_ = new BSTNode_();
	init_bst_node_(root_, root_id, clustering_heuristic_->get_cluster_radius(root_id), graph_ids);

	// Recursively build the BST and mark the BST as initialized.
	if (root_->cluster.size() > max_cluster_size) {
		ProgressBar progress(graph_ids.size());
		std::cout << "Initializing BST: " << progress << std::flush;
		split_(max_cluster_size, root_, pos_next_focal_graph_id, progress);
		std::cout << "\n";
	}

	// Set the initialization time.
	auto end = std::chrono::high_resolution_clock::now();
	init_time_ = start - end;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
init_bst_node_(BSTNode_ * bst_node, GEDGraph::GraphID focal_graph_id, double radius, const std::vector<GEDGraph::GraphID> & cluster) {
	bst_node->focal_graph_id = focal_graph_id;
	bst_node->radius = radius;
	bst_node->cluster = cluster;
	bst_node->child_1 = nullptr;
	bst_node->child_2 = nullptr;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
std::size_t
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
split_(std::size_t max_cluster_size, BSTNode_ * node, std::size_t pos_next_focal_graph_id, ProgressBar & progress) {

	// Get the IDs of the children's focal graphs.
	GEDGraph::GraphID child_id_1{focal_graph_ids_.at(pos_next_focal_graph_id++)};
	GEDGraph::GraphID child_id_2{focal_graph_ids_.at(pos_next_focal_graph_id++)};

	// Run the graph clustering heuristic and obtain the clustering.
	clustering_heuristic_->run(node->cluster, {child_id_1, child_id_2});
	std::map<GEDGraph::GraphID, std::vector<GEDGraph::GraphID>> clustering;
	clustering_heuristic_->get_clustering(clustering);

	// Initialize the children.
	node->child_1 = new BSTNode_();
	init_bst_node_(node->child_1, child_id_1, clustering_heuristic_->get_cluster_radius(child_id_1), clustering.at(child_id_1));
	node->child_2 = new BSTNode_();
	init_bst_node_(node->child_2, child_id_2, clustering_heuristic_->get_cluster_radius(child_id_2), clustering.at(child_id_2));

	// Print information.
	progress.increment();
	std::cout << "Initializing BST: " << progress << std::flush;

	// Recursively call split_() on the children.
	if (node->child_1->cluster.size() > max_cluster_size) {
		pos_next_focal_graph_id = split_(max_cluster_size, node->child_1, pos_next_focal_graph_id, progress);
	}
	if (node->child_2->cluster.size() > max_cluster_size) {
		pos_next_focal_graph_id = split_(max_cluster_size, node->child_2, pos_next_focal_graph_id, progress);
	}

	// Return the position of the next focal graph ID.
	return pos_next_focal_graph_id;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
process_range_query(GEDGraph::GraphID query_graph_id, double threshold) {

	std::cout << "\n===========================================================\n";
	std::cout << "Processing range query.\n";
	std::cout << "-----------------------------------------------------------\n";

	// Ensure that the BST has been initialized.
	if (root_ == nullptr) {
		throw Error("The GraphBST has not been initialized. Call init() or load() before calling process_range_query().");
	}

	// Start recording the query time.
	auto start = std::chrono::high_resolution_clock::now();

	// Clear the result variables.
	verified_graph_ids_.clear();
	filtered_graph_ids_.clear();
	undecided_graph_ids_.clear();
	num_ged_evals_ = 0;

	// Recursively process the query.
	check_(query_graph_id, threshold, root_);

	// Set the query time.
	auto end = std::chrono::high_resolution_clock::now();
	query_time_ = end - start;

	std::cout << "Number of GED evaluations: " << num_ged_evals_ << "\n";
	std::cout << "Number of verified data graphs: " << verified_graph_ids_.size() << "\n";
	std::cout << "Number of filtered data graphs: " << filtered_graph_ids_.size() << "\n";
	std::cout << "Number of undecided data graphs: " << undecided_graph_ids_.size() << "\n";
	std::cout << "Initialization time: " << get_init_time() << "\n";
	std::cout << "Query time: " << get_query_time() << "\n";
	std::cout << "===========================================================\n";
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
check_(GEDGraph::GraphID query_graph_id, double threshold, const BSTNode_ * node) {

	// Update number of GED evaluations.
	num_ged_evals_++;

	// Compute lower bound.
	ged_env_->set_method(lower_bound_method_, lower_bound_options_);
	ged_env_->run_method(query_graph_id, node->focal_graph_id);

	// Use lower bound for filtering.
	if (ged_env_->get_lower_bound(query_graph_id, node->focal_graph_id) > threshold + node->radius) {
		for (GEDGraph::GraphID graph_id : node->cluster) {
			filtered_graph_ids_.emplace_back(graph_id);
		}
		return;
	}

	// Compute upper bound if the upper bound method is different from the lower bound method.
	if (lower_bound_method_ != upper_bound_method_ or lower_bound_options_ != upper_bound_options_) {
		ged_env_->set_method(upper_bound_method_, upper_bound_options_);
		ged_env_->run_method(query_graph_id, node->focal_graph_id);
	}

	// Use upper bound for verifying.
	if (ged_env_->get_upper_bound(query_graph_id, node->focal_graph_id) <= threshold - node->radius) {
		for (GEDGraph::GraphID graph_id : node->cluster) {
			verified_graph_ids_.emplace_back(graph_id);
		}
		return;
	}

	// The current node is a leaf that has neither been filtered nor verified.
	if (node->child_1 == nullptr) {
		num_ged_evals_ += node->cluster.size();
		for (GEDGraph::GraphID graph_id : node->cluster) {
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
	check_(query_graph_id, threshold, node->child_1);
	check_(query_graph_id, threshold, node->child_2);

}

template<>
void
GraphBST<GXLNodeID, GXLLabel, GXLLabel>::
serialize_(const BSTNode_ * node, const std::string & node_id, const std::string & focal_graph_dir, std::map<std::string, std::string> & config) const {
	// Write "NULL" and return if the current node equals nullptr.
	if (node == nullptr) {
		config[node_id] = "NULL";
		return;
	}

	// Save the focal graph of the current node.
	ged_env_->save_as_gxl_graph(node->focal_graph_id, focal_graph_dir + "/" + ged_env_->get_graph_name(node->focal_graph_id));

	// Save the radius and the cluster in the configuration map.
	config[node_id] = ged_env_->get_graph_name(node->focal_graph_id) + "," + std::to_string(node->radius);
	for (GEDGraph::GraphID graph_id : node->cluster) {
		config[node_id] += "," + ged_env_->get_graph_name(graph_id);
	}

	// Continue with the children.
	serialize_(node->child_1, node_id + "1", focal_graph_dir, config);
	serialize_(node->child_2, node_id + "2", focal_graph_dir, config);
}

template<>
void
GraphBST<GXLNodeID, GXLLabel, GXLLabel>::
save(const std::string & bst_file_name, const std::string & focal_graph_dir) const {
	// Serialize the tree.
	std::map<std::string, std::string> config;
	serialize_(root_, "0", focal_graph_dir, config);

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
	util::tokenize(config.at(node_id), ',', node_info);

	// Load the cluster and load the data graphs into the environment, if the current node is the root.
	std::vector<GEDGraph::GraphID> cluster;
	for (std::size_t pos{2}; pos < node_info.size(); pos++) {
		if (node_id == "0") {
			GEDGraph::GraphID graph_id{ged_env_->load_gxl_graph(data_graph_dir, node_info.at(pos), node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, undefined(), "no_class")};
			graph_name_to_id[node_info.at(pos)] = graph_id;
			cluster.push_back(graph_id);
		}
		else {
			cluster.push_back(graph_name_to_id.at(node_info.at(pos)));
		}
	}

	// Load the focal graph into the environment.
	GEDGraph::GraphID focal_graph_id{ged_env_->load_gxl_graph(focal_graph_dir, node_info.at(0), node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, undefined(), "no_class")};

	// Get the radius.
	double radius{std::stod(node_info.at(1))};

	// Initialize the current node.
	node = new BSTNode_();
	init_bst_node_(node, focal_graph_id, radius, cluster);

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
	de_serialize_("0", config, data_graph_dir, focal_graph_dir, node_type, edge_type, irrelevant_node_attributes, irrelevant_edge_attributes, root_, graph_name_to_id);
	ged_env_->init(ged_env_->get_init_type());

	// Set the initialization time.
	auto end = std::chrono::high_resolution_clock::now();
	init_time_ = end - start;
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
	return num_ged_evals_;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
BSTNode_::
BSTNode_() :
focal_graph_id{undefined()},
radius{std::numeric_limits<double>::infinity()},
cluster(),
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
