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
 * @file graph_clustering_heuristic.hpp
 * @brief ged::GraphClusteringHeuristic class declaration.
 */

#ifndef MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTICS_HPP_
#define MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTICS_HPP_

#include "median_graph_estimator.hpp"

namespace ged {

/*!
 * @brief Class for clustering a collection of graphs.
 *
 * @details Implements the algorithm suggested in:
 * - N. Boria, S. Bougleux, B. Ga&uuml;z&egrave;re, D. B. Blumenthal, and L. Brun:
 *   &ldquo;Scalable generalized graph median estimation and applications in clustering, classification, and indexing&rdquo;
 *   submitted to VLDB J.,
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--focal-graphs MEDIANS\|MEDOIDS</tt> | use medians or medoids as the focal graphs of the clusters | @p MEDIANS | n.a. |
 * | <tt>\--init-type CLUSTERS\|K-MEANS++</tt> | approach used for generating initial clusters | @p K-MEANS++ | if @p K-MEANS++, well distributed graphs are used as the first focal graphs |
 * | <tt>\--random-inits @<convertible to int greater 0@></tt> | number of randomly constructed initial clusterings | @p 50 | if @p 1, the option @p \--minimize has no effect |
 * | <tt>\--randomness REAL\|PSEUDO</tt> | use real randomness or pseudo randomness | @p REAL | if @p REAL, the option @p \--seed has no effect |
 * | <tt>\--seed @<convertible to int greater equal 0@></tt> | seed for generating pseudo random numbers | @p 0 | n.a. |
 * | <tt>\--max-itrs @<convertible to int@></tt> | maximal number of iterations in Lloyd's algorithm | @p 100 | if negative, no maximal number of iterations is enforced |
 * | <tt>\--time-limit @<convertible to double@></tt> | time limit in seconds for main block gradient descent | @p 0 | if less or equal @p 0, no time limit is enforced |
 * | <tt>\--epsilon @<convertible to double greater 0@></tt> | convergence threshold used everywhere | @p 0.0001 | n.a. |
 * | <tt>\--stdout 0\|1\|2</tt> | print runtime information to standard output stream | @p 2 | @p 0: no output; @p 1: output only before termination; @p 2: output also during optimization |
 */
template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
class GraphClusteringHeuristic {

public:

	/*!
	 * @brief Constructor.
	 * @param[in] ged_env Pointer to initialized environment. The edit costs must be set by the user.
	 * @param[in] mge Pointer to median graph estimator constructed on top of the environment @p ged_env.
	 * You can set this argument to nullptr if you always use the clustering heuristic with one of the
	 * options "--max-itrs 0" (random clustering), "--clustering-method MEDOID" or "--clustering-method CENTER".
	 */
	GraphClusteringHeuristic(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge);

	/*!
	 * @brief Sets the options of the clustering heuristic.
	 * @param[in] options String that specifies with which options to run the clustering heuristic.
	 */
	void set_options(const std::string & options);

	/*!
	 * @brief Selects method to be used for computing (upper bounds for) GED.
	 * @param[in] ged_method The selected method. Default: ged::Options::GEDMethod::BRANCH_UNIFORM.
	 * @param[in] ged_options The options for the selected method. Default: "".
	 */
	void set_ged_method(Options::GEDMethod ged_method, const std::string ged_options);

	/*!
	 * @brief Runs the graph clustering algorithm.
	 * @param[in] graph_ids Vector that contains the IDs of the graphs that should be clustered.
	 * @param[in] focal_graph_ids Vector that contains the IDs of the cluster's focal graph.
	 * @note The number of clusters equals the lenght of @p focal_graph_ids.
	 */
	void run(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids);

	/*!
	 * @brief Returns the runtime.
	 * @return Runtime in seconds.
	 */
	double get_runtime() const;

	/*!
	 * @brief Get the clustering.
	 * @param[out] clustering A map that contains pairs of the form <tt>(focal_graph_id, cluster)</tt>,
	 * where @p cluster contains the IDs of the graphs that are contained in the cluster associated to the focal graph with ID @p focal_graph_id.
	 */
	void get_clustering(std::map<GEDGraph::GraphID, std::vector<GEDGraph::GraphID>> & clustering) const;

	/*!
	 * @brief Return overall sum of distances.
	 * @return Overall sum of distances.
	 */
	double get_sum_of_distances() const;

	/*!
	 * @brief Returns cluster radius.
	 * @param[in] focal_graph_id ID of the focal graph whose cluster radius should be returned.
	 * @return Radius of the cluster associated to the focal graph with ID @p focal_graph_id.
	 */
	double get_cluster_radius(GEDGraph::GraphID focal_graph_id) const;

	/*!
	 * @brief Returns sum of distances within cluster.
	 * @param[in] focal_graph_id ID of the focal graph for which the sum of distances of the cluster should be returned.
	 * @return Sum of distances of the cluster associated to the focal graph with ID @p focal_graph_id.
	 */
	double get_cluster_sum_of_distances(GEDGraph::GraphID focal_graph_id) const;

	/*!
	 * @brief Returns the ID of the assigned focal graph (a.k.a. cluster ID).
	 * @param[in] graph_id ID of the graph for which the ID of the assigned focal graph should be returned.
	 * @return ID of the assigned focal graph.
	 */
	GEDGraph::GraphID get_assigned_focal_graph_id(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Returns the distance from the assigned focal graph.
	 * @param[in] graph_id ID of the graph whose distance from its assigned focal graph should be returned.
	 * @return Distance from the assigned focal graph to the graph with ID @p graph_id.
	 */
	double get_distance_from_assigned_focal_graph(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Returns the node map from the assigned focal graph.
	 * @param[in] graph_id ID of the graph whose node map from its assigned focal graph should be returned.
	 * @return Node map from the assigned focal graph to the graph with ID @p graph_id.
	 */
	const NodeMap & get_node_map_from_assigned_focal_graph(GEDGraph::GraphID graph_id) const;

	/*!
	 * @brief Saves the computed focal graphs as GXL graphs and creates a GraphCollection file that lists all of them.
	 * @param[in] collection_file_name The name of the collection file.
	 * @param[in] focal_graph_file_names A map that contains the names of the GXL files for the focal graphs.
	 * Must contain a pair <tt>(focal_graph_id, focal_graph_file_name)</tt> for each @p focal_graph_id passed to run().
	 * @param[in] focal_graph_classes A vector that contains the classes of the focal graphs. If left empty, no classes are specified in the collection file.
	 * Otherwise, it must have the same length as the argument @p focal_graph_ids passed to run().
	 */
	void save_focal_graphs(const std::string & collection_file_name, const std::map<GEDGraph::GraphID, std::string> & focal_graph_file_names, const std::map<GEDGraph::GraphID, std::string> & focal_graph_classes = {}) const;

	/*!
	 * @brief Computes the normalized mutual information between the computed clustering and a ground truth clustering.
	 * @param[in] ground_truth_clustering The ground truth clustering.
	 * @return Normalized mutual information between the two clusterings, i.e., a score between 0 and 1 that equals 1 just in case the two clusterings are identical.
	 */
	double get_normalized_mutual_information(const std::vector<std::vector<GEDGraph::GraphID>> & ground_truth_clustering) const;

	/*!
	 * @brief Computes the adjusted Rand index between the computed clustering and a ground truth clustering.
	 * @param[in] ground_truth_clustering The ground truth clustering.
	 * @return Adjusted Rand index between the two clusterings, i.e., a score between 0 and 1 that equals 1 just in case the two clusterings are identical.
	 */
	double get_adjusted_rand_index(const std::vector<std::vector<GEDGraph::GraphID>> & ground_truth_clustering) const;

private:

	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env_;

	MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge_;

	Options::GEDMethod ged_method_;

	std::string ged_options_;

	std::string clustering_method_;

	std::string init_type_;

	bool use_real_randomness_;

	std::size_t num_random_inits_;

	std::size_t seed_;

	double time_limit_in_sec_;

	double epsilon_;

	std::size_t max_itrs_;

	std::size_t print_to_stdout_;

	std::map<GEDGraph::GraphID, GEDGraph::GraphID> assigned_focal_graph_ids_;

	std::map<GEDGraph::GraphID, NodeMap> node_maps_from_assigned_focal_graphs_;

	std::map<GEDGraph::GraphID, double> cluster_sums_of_distances_;

	std::map<GEDGraph::GraphID, double> cluster_radii_;

	double sum_of_distances_;

	Seconds runtime_;

	std::vector<std::size_t> itrs_;

	void set_default_options_();

	void initialize_focal_graphs_and_clusters_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
			std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs);

	void compute_initial_focal_graphs_k_means_plus_plus_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
			std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs);

	void compute_initial_focal_graphs_cluster_sampling_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
			std::mt19937 & urng, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs);

	bool termination_criterion_met_(bool converged, const Timer & timer, std::size_t itr) const;

	bool update_clusters_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids);

	void update_focal_graphs_(const std::vector<GEDGraph::GraphID> & focal_graph_ids, std::map<GEDGraph::GraphID, ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel>> & focal_graphs);

	double compute_medoid_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID medoid_id,
			ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & medoid, std::map<GEDGraph::GraphID, NodeMap> & node_maps_from_medoid) const;

	double compute_median_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID median_id,
			ExchangeGraph<UserNodeID, UserNodeLabel, UserEdgeLabel> & median, std::map<GEDGraph::GraphID, NodeMap> & node_maps_from_median) const;

	void compute_cluster_radii_(const std::vector<GEDGraph::GraphID> & focal_graph_ids);

	double compute_intersection_size_(const std::vector<GEDGraph::GraphID> & sorted_cluster_1, const std::vector<GEDGraph::GraphID> & sorted_cluster_2) const;

};

}

#include "graph_clustering_heuristic.ipp"

#endif /* MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTICS_HPP_ */
