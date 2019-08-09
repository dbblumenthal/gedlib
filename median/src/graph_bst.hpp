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
 * @file graph_bst.hpp
 * @brief ged::GraphBST class declaration.
 */

#ifndef MEDIAN_SRC_GRAPH_BST_HPP_
#define MEDIAN_SRC_GRAPH_BST_HPP_

#include "median_graph_estimator.hpp"

namespace ged {

/*!
 * @brief Class for clustering a collection of graphs.
 *
 * @details Implements monotonic graph bisector trees as suggested in:
 * - N. Boria, S. Bougleux, B. Ga&uuml;z&egrave;re, D. B. Blumenthal, and L. Brun:
 *   &ldquo;Scalable generalized graph median estimation and applications in clustering, classification, and indexing&rdquo;
 *   submitted to VLDB J.,
 *
 *   Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--focal-graphs MEDIANS\|MEDOIDS\|CENTERS</tt> | use medians, medoids, or centers as the focal graphs of the tree's inner nodes | @p MEDIANS | n.a. |
 * | <tt>\--max-cluster-size @<convertible to int greater 0@></tt> | maximal size of the clusters of the tree's inner nodes | @p 10 | n.a. |
 * | <tt>\--cutoff @<convertible to double greater 0 and smaller equal 1@></tt> | employed to determine which graphs are used to compute the new focal graphs when splitting inner nodes | @p 0.5 | n.a. |
 * | <tt>\--stdout 0\|1\|2</tt> | print runtime information to standard output stream | @p 2 | @p 0: no output; @p 1: output only before termination; @p 2: output also during optimization |
 */
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
class GraphBST {

public:

	/*!
	 * @brief Constructor.
	 * @param[in] ged_env Pointer to initialized environment. The edit costs must be set by the user.
	 * @param[in] clustering_heuristic Pointer to graph clustering heuristic constructed on top of the environment @p ged_env.
	 */
	GraphBST(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge);

	/*!
	 * @brief Destructor.
	 */
	~GraphBST();

	/*!
	 * @brief Sets the options of the estimator.
	 * @param[in] options String that specifies with which options to run the estimator.
	 */
	void set_options(const std::string & options);

	/*!
	 * @brief Selects the lower bound method to be used for processing range queries.
	 * @param[in] lower_bound_method The selected method. Default: ged::Options::GEDMethod::BRANCH_FAST.
	 * @param[in] lower_bound_options The options for the selected main method. Default: "".
	 */
	void set_lower_bound_method(Options::GEDMethod lower_bound_method, const std::string & lower_bound_options = "");

	/*!
	 * @brief Selects the upper bound method to be used for processing range queries.
	 * @param[in] upper_bound_method The selected method. Default: ged::Options::GEDMethod::IPFP.
	 * @param[in] upper_bound_options The options for the selected main method. Default: "".
	 */
	void set_upper_bound_method(Options::GEDMethod upper_bound_method, const std::string & upper_bound_options = "");

	/*!
	 * @brief Initializes the bisector tree.
	 * @param[in] graph_ids Vector that contains the IDs of the data graphs that should be indexed.
	 * @param[in] focal_graph_ids Vector that contains the IDs of the focal graphs used at the inner nodes of the tree.
	 * Must be of size at least <tt>graph_ids.size()</tt>.
	 */
	void init(std::vector<GEDGraph::GraphID> graph_ids, std::vector<GEDGraph::GraphID> focal_graph_ids);

	/*!
	 * @brief Saves the bisector tree.
	 * @param[in] bst_file_name The name of the configuration file where the tree should be stored.
	 * @param[in] focal_graph_dir The directory where the focal graphs of the tree's inner nodes should be stored. If empty string, only the tree is saved.
	 * @param[in] max_cluster_size Clusters of size at most @p max_cluster_size are saved as leafs. If @p 0, the internal tree is saved as it is.
	 * @param[in] only_tree
	 */
	void save(const std::string & bst_file_name, const std::string & focal_graph_dir, std::size_t max_cluster_size = 0) const;

	/*!
	 * @brief Load the bisector tree from a configuration file.
	 * @param[in] bst_file_name The name of the configuration file.
	 * @param[in] data_graph_dir The directory of the data graphs.
	 * @param[in] focal_graph_dir The directory of the focal graphs of the tree's inner nodes.
	 * @param[in] node_type Select if nodes are labeled or unlabeled.
	 * @param[in] edge_type Select if edges are labeled or unlabeled.
	 * @param[in] irrelevant_node_attributes Set of node attributes that are irrelevant for the selected edit costs.
	 * @param[in] irrelevant_edge_attributes Set of edge attributes that are irrelevant for the selected edit costs.
	 */
	void load(const std::string & bst_file_name, const std::string & data_graph_dir, const std::string & focal_graph_dir,
			Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
			const std::unordered_set<std::string> & irrelevant_node_attributes = {}, const std::unordered_set<std::string> & irrelevant_edge_attributes = {});

	/*!
	 * @brief Process a GED range query.
	 * @param[in] query_graph_id The ID of the query graph.
	 * @param[in] threshold The range threshold.
	 * @warning Throws an error if neither init() nor load() has been called before.
	 */
	void process_range_query(GEDGraph::GraphID query_graph_id, double threshold);

	/*!
	 * @brief Returns the IDs of the graphs that have been verified.
	 * @return IDs of the graphs that have been verified.
	 */
	const std::vector<GEDGraph::GraphID> & get_verified_graph_ids() const;

	/*!
	 * @brief Returns the IDs of the graphs that have been filtered.
	 * @return IDs of the graphs that have been filtered.
	 */
	const std::vector<GEDGraph::GraphID> & get_filtered_graph_ids() const;

	/*!
	 * @brief Returns the IDs of the graphs that have been neither filtered nor verified.
	 * @return IDs of the graphs that have been neither filtered nor verified.
	 */
	const std::vector<GEDGraph::GraphID> & get_undecided_graph_ids() const;

	/*!
	 * @brief Returns the initialization time.
	 * @return Runtime of last call to init() or load() in seconds.
	 */
	double get_init_time() const;

	/*!
	 * @brief Returns the query time.
	 * @return Runtime of the last call to process_range_query().
	 */
	double get_query_time() const;

	/*!
	 * @brief Returns total the number of GED evaluations.
	 * @return The total number of GED evaluations.
	 */
	std::size_t get_num_ged_evals() const;

	/*!
	 * @brief Returns the number of GED evaluations between the query graph and the focal graphs stored at the inner nodes of the tree.
	 * @return The number of GED evaluations between the query graph and the focal graphs stored at the inner nodes of the tree.
	 */
	std::size_t get_num_ged_evals_focal_graphs() const;

	/*!
	 * @brief Returns the number of GED evaluations between the query graph and the graphs stored in the leafs of the tree.
	 * @return The number of GED evaluations between the query graph and and the graphs stored in the leafs of the tree.
	 */
	std::size_t get_num_ged_evals_data_graphs() const;

private:

	struct BSTNode_ {

		BSTNode_();

		~BSTNode_();

		GEDGraph::GraphID focal_graph_id;

		double radius;

		std::vector<GEDGraph::GraphID> cluster;

		std::map<GEDGraph::GraphID, double> distances_from_focal_graph;

		GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::BSTNode_ * child_1;

		GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::BSTNode_ * child_2;
	};

	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env_;

	MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge_;

	Options::GEDMethod lower_bound_method_;

	std::string lower_bound_options_;

	Options::GEDMethod upper_bound_method_;

	std::string upper_bound_options_;

	std::size_t max_cluster_size_;

	std::string focal_graphs_;

	double cutoff_;

	std::size_t print_to_stdout_;

	BSTNode_ * root_;

	std::vector<GEDGraph::GraphID> focal_graph_ids_;

	std::vector<GEDGraph::GraphID> verified_graph_ids_;

	std::vector<GEDGraph::GraphID> filtered_graph_ids_;

	std::vector<GEDGraph::GraphID> undecided_graph_ids_;

	Seconds init_time_;

	Seconds query_time_;

	std::size_t num_ged_evals_focal_graphs_;

	std::size_t num_ged_evals_data_graphs_;

	void set_default_options_();

	void init_bst_node_(BSTNode_ * bst_node, GEDGraph::GraphID focal_graph_id, double radius, const std::vector<GEDGraph::GraphID> & cluster, const std::map<GEDGraph::GraphID, double> & distances_from_focal_graph);

	std::size_t split_(BSTNode_ * node, std::size_t pos_next_focal_graph_id, ProgressBar & progress);

	void compute_new_focal_graph_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID new_focal_graph_id, std::map<GEDGraph::GraphID, double> & distances_from_new_focal_graph);

	void check_(GEDGraph::GraphID query_graph_id, double threshold, const BSTNode_ * node, double lower_bound, double upper_bound);

	void serialize_(const BSTNode_ * node, const std::string & node_id, std::size_t max_cluster_size, const std::string & focal_graph_dir, std::map<std::string, std::string> & config) const;

	void de_serialize_(const std::string & node_id, const std::map<std::string, std::string> & config, const std::string & data_graph_dir, const std::string & focal_graph_dir,
			Options::GXLNodeEdgeType node_type, Options::GXLNodeEdgeType edge_type,
			const std::unordered_set<std::string> & irrelevant_node_attributes, const std::unordered_set<std::string> & irrelevant_edge_attributes,
			BSTNode_ * & node, std::map<std::string, GEDGraph::GraphID> & graph_name_to_id);

};

}

#include "graph_bst.ipp"

#endif /* MEDIAN_SRC_GRAPH_BST_HPP_ */
