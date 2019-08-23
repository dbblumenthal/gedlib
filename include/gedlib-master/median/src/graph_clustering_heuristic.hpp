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
 * | <tt>\--clustering-method K-MEDIANS\|K-MEDOIDS</tt> | method for computing the focal graphs of the clusters | @p MEDIAN | n.a. |
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

	void run(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids);

	void get_runtime() const;

	void get_clustering(std::vector<std::vector<GEDGraph::GraphID>> & clustering) const;

	double get_sum_of_distances() const;

	double get_cluster_radius(GEDGraph::GraphID median_id) const;

	double get_cluster_sum_of_distances(GEDGraph::GraphID median_id) const;

	GEDGraph::GraphID get_assigned_focal_graph_id(GEDGraph::GraphID graph_id) const;

	double get_distance_from_assigned_focal_graph(GEDGraph::GraphID graph_id) const;

	const NodeMap & get_node_map_from_assigned_focal_graph(GEDGraph::GraphID graph_id) const;

	double get_clustering_distance(const std::vector<std::vector<GEDGraph::GraphID>> & ground_truth_clustering) const;

	void save_focal_graphs(const std::string & collection_file_name, const std::vector<std::string> & focal_graph_file_names, const std::vector<std::string> & focal_graph_classes = {}) const;

private:

	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env_;

	MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge_;

	Options::GEDMethod ged_method_;

	std::string ged_options_;

	std::string clustering_method_;

	std::string init_type_;

	bool use_real_randomness_;

	std::size_t seed_;

	std::size_t num_random_inits_;

	double time_limit_in_sec_;

	double epsilon_;

	std::size_t max_itrs_;

	std::size_t print_to_stdout_;

	std::map<GEDGraph::GraphID, GEDGraph::GraphID> assigned_focal_graph_ids_;

	std::map<GEDGraph::GraphID, NodeMap> node_maps_from_assigned_focal_graphs_;

	std::map<GEDGraph::GraphID, double> cluster_sums_of_distances_;

	std::map<GEDGraph::GraphID, double> cluster_radii_;

	double sum_of_distances_;

	std::vector<std::size_t> itrs_;

	Seconds runtime_;

	void set_default_options_();

	void compute_initial_clusters_(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & focal_graph_ids,
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
