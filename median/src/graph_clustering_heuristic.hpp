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
 * | <tt>\--random-inits @<convertible to int greater 0@></tt> | number of randomly constructed initial clusterings | @p 10 | n.a. |
 * | <tt>\--randomness REAL\|PSEUDO</tt> | use real randomness or pseudo randomness | @p REAL | if @p REAL, the option @p \--seed has no effect |
 * | <tt>\--seed @<convertible to int greater equal 0@></tt> | seed for generating pseudo random numbers | @p 0 | n.a. |
 * | <tt>\--refine TRUE\|FALSE</tt> | improve node maps and sums of distances for converged clusters | @p TRUE | n.a. |
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
	 * options "--max-itrs 0" (random clustering) or "--focal-graphs MEDOIDS".
	 */
	GraphClusteringHeuristic(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge);

	/*!
	 * @brief Sets the options of the clustering heuristic.
	 * @param[in] options String that specifies with which options to run the clustering heuristic.
	 */
	void set_options(const std::string & options);

	/*!
	 * @brief Selects the method used during the computation of the final clusters.
	 * @param[in] main_method The selected method. Default: ged::Options::GEDMethod::BRANCH_FAST.
	 * @param[in] main_options The options for the selected main method. Default: "".
	 */
	void set_main_method(Options::GEDMethod main_method, const std::string & main_options = "");

	/*!
	 * @brief Selects method to be used for improving the sums of distances of the converged clusters.
	 * @param[in] refine_method The selected method. Default: ged::Options::GEDMethod::IPFP.
	 * @param[in] refine_options The options for the selected method. Default: "".
	 * @note Has no effect if "--refine FALSE" is passed to set_options().
	 */
	void set_refine_method(Options::GEDMethod refine_method, const std::string & refine_options = "");

	/*!
	 * @brief Runs the graph clustering algorithm.
	 * @param[in] graph_ids Vector that contains the IDs of the graphs that should be clustered.
	 * @param[in] focal_graph_ids Vector that contains the IDs of the clusters' focal graphs.
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
	 * @brief Returns number of iterations.
	 * @return A vector that contains the number of iterations for each initial median for the last call to run().
	 */
	const std::vector<std::size_t> & get_num_itrs() const;

	/*!
	 * @brief Saves the computed focal graphs as GXL graphs and creates a GraphCollection file that lists all of them.
	 * @param[in] collection_file_name The name of the collection file.
	 * @param[in] focal_graph_dir The directory where the focal graphs should be stored as GXL files.
	 */
	void save(const std::string & collection_file_name, const std::string & focal_graph_dir) const;

	/*!
	 * @brief Computes the adjusted Rand index between the computed clustering and a ground truth clustering.
	 * @param[in] ground_truth_clustering The ground truth clustering.
	 * @return Adjusted Rand index between the two clusterings, i.e., a score between -1 and 1 that equals 1 just in case the two clusterings are identical.
	 */
	double get_adjusted_rand_index(const std::vector<std::vector<GEDGraph::GraphID>> & ground_truth_clustering) const;

	/*!
	 * @brief Computes the Gini coefficient of the computed clustering.
	 * @return Gini coefficient of the computed clustering, i.e., a score between 0 and 1 that equals 0 just in case the clustering is completely balanced and 1 just in case one cluster contains all data graphs.
	 */
	double get_gini_coefficient() const;
    
    /*!
     * @brief Computes the mean silhouette coefficient of all clustered graphs.
     * @return The mean silhouette score of all clustered graphs.
     * The silhouette a of a clustered graph is a score between -1 and 1 that is close to 1 if the mean distance from the graph to the other graphs contained in its own cluster is much smaller than the mean distance
     * to the graphs in the closest other cluster.
     */
    double get_silhouette_score() const;

	/*!
	 * @brief Returns pointer to the environment employed by the clustering heuristic.
	 * @return Pointer to the environment employed by the clustering heuristic.
	 */
	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * get_ged_env();

private:

	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env_;

	MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge_;

	Options::GEDMethod main_method_;

	std::string main_options_;

	Options::GEDMethod refine_method_;

	std::string refine_options_;

	std::string focal_graphs_;

	std::string init_type_;

	bool use_real_randomness_;

	std::size_t num_random_inits_;

	std::size_t seed_;

	bool refine_;

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

	bool update_clusters_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids, bool refine);

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
