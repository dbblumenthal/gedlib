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
query_time_() {
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
	if (not root_ == nullptr) {
		delete root_;
		root_ = nullptr;
	}
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
init(std::vector<GEDGraph::GraphID> graph_ids, std::size_t max_cluster_size) {

	// Start recording the initialization time.
	auto start = std::chrono::high_resolution_clock::now();

	// Unitialize the BST if it has already been initialized.
	if (not root_ == nullptr) {
		delete root_;
		root_ = nullptr;
	}

	// Allocate space for the focal graphs and re-initialize the environment.
	focal_graph_ids_.clear();
	for (std::size_t counter{0}; 0 < 2 * graph_ids.size(); counter++) {
		focal_graph_ids_.emplace_back(ged_env_->add_graph("focal_graph_" + std::to_string(counter) + ".gxl", "no_class"));
	}
	ged_env_->init(ged_env_->get_init_type());

	// Initialize the root of the BST.
	std::size_t pos_next_focal_graph_id{0};
	GEDGraph::GraphID root_id{focal_graph_ids_.at(pos_next_focal_graph_id++)};
	clustering_heuristic_->run(graph_ids, {root_id});
	root_ = new BSTNode_();
	init_bst_node_(root_, root_id, clustering_heuristic_->get_cluster_radius(root_id), graph_ids);

	// Recursively build the BST and mark the BST as initialized.
	if (root_->cluster.size() > max_cluster_size) {
		split_(max_cluster_size, root_, pos_next_focal_graph_id);
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
split_(std::size_t max_cluster_size, const BSTNode_ * node, std::size_t pos_next_focal_graph_id) {

	// Get the IDs of the children's focal graphs.
	GEDGraph::GraphID child_id_1{focal_graph_ids_.at(pos_next_focal_graph_id++)};
	GEDGraph::GraphID child_id_2{focal_graph_ids_.at(pos_next_focal_graph_id++)};

	// Run the graph clustering heuristic and obtain the clustering.
	clustering_heuristic_->run(node->cluster(), {child_id_1, child_id_2});
	std::map<GEDGraph::GraphID, std::vector<GEDGraph::GraphID>> clustering;
	clustering_heuristic_->get_clustering(clustering);

	// Initialize the children.
	node->child_1 = new BSTNode_();
	init_bst_node_(node->child_1, child_id_1, clustering_heuristic_->get_cluster_radius(child_id_1), clustering.at(child_id_1));
	node->child_2 = new BSTNode_();
	init_bst_node_(node->child_2, child_id_2, clustering_heuristic_->get_cluster_radius(child_id_2), clustering.at(child_id_2));

	// Recursively call split_() on the children.
	if (node->child_1->cluster.size() > max_cluster_size) {
		pos_next_focal_graph_id = split_(max_cluster_size, node->child_1, pos_next_focal_graph_id);
	}
	if (node->child_2->cluster.size() > max_cluster_size) {
		pos_next_focal_graph_id = split_(max_cluster_size, node->child_2, pos_next_focal_graph_id);
	}

	// Return the position of the next focal graph ID.
	return pos_next_focal_graph_id;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
process_range_query(GEDGraph::GraphID query_graph_id, double threshold) {

	// Ensure that the BST has been initialized.
	if (root_ == nullptr) {
		throw Error("The GraphBST has not been initialized. Call init() or load() before calling process_range_query().");
	}

	// Start recording the query time.
	auto start = std::chrono::high_resolution_clock::now();

	// Set the threshold and clear the result sets.
	verified_graph_ids_.clear();
	filtered_graph_ids_.clear();
	undecided_graph_ids_.clear();

	// Recursively process the query.
	check_(query_graph_id, threshold, root_);

	// Set the query time.
	auto end = std::chrono::high_resolution_clock::now();
	query_time_ = start - end;
}

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
void
GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::
check_(GEDGraph::GraphID query_graph_id, double threshold, const BSTNode_ * node) {

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


}


#endif /* MEDIAN_SRC_GRAPH_BST_IPP_ */
