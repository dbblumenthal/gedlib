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

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
class GraphClusteringHeuristic {

public:

	/*!
	 * @brief Constructor.
	 * @param[in] ged_env Pointer to initialized environment. The edit costs must be set by the user.
	 * @param[in] mge Pointer to median graph estimator constructed on top of the environment @p ged_env.
	 * You can set this argument to nullptr if you always use the clustering heuristic with one of the options "--clustering-method RANDOM" or "--clustering-method MEDOID".
	 */
	GraphClusteringHeuristic(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge);

	void set_options(const std::string & options);

	void set_medoid_method(Options::GEDMethod medoid_method, const std::string medoid_options);

	void run(const std::vector<GEDGraph::GraphID> & graph_ids, const std::vector<GEDGraph::GraphID> & median_ids);

	void get_runtime() const;

	void get_clustering(std::vector<std::vector<GEDGraph::GraphID>> & clustering) const;

	double get_sum_of_distances() const;

	double get_cluster_radius(GEDGraph::GraphID median_id) const;

	double get_cluster_sum_of_distances(GEDGraph::GraphID median_id) const;

	GEDGraph::GraphID get_assigned_median_id(GEDGraph::GraphID graph_id) const;

	double get_distance_from_assigned_median(GEDGraph::GraphID graph_id) const;

	const NodeMap & get_node_map_from_assigned_median(GEDGraph::GraphID graph_id) const;

	double get_clustering_distance(const std::vector<std::vector<GEDGraph::GraphID>> & ground_truth_clustering) const;

	void save_medians(const std::string & collection_file_name, const std::vector<std::string> & median_file_names, const std::vector<std::string> & median_classes = {}) const;

private:

	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env_;

	MedianGraphEstimator<UserNodeID, UserNodeLabel, UserEdgeLabel> * mge_;

	Options::GEDMethod medoid_method_;

	std::string medoid_options_;

	std::string clustering_method_;

	std::string init_type_;

	bool use_real_randomness_;

	std::size_t seed_;

	std::size_t num_inits_;

	std::size_t max_itrs_;

	double time_limit_in_sec_;

	std::size_t print_to_stdout_;

	std::map<GEDGraph::GraphID, GEDGraph::GraphID> assigned_medians_;

	std::map<GEDGraph::GraphID, NodeMap> node_maps_from_assigned_medians_;

	std::map<GEDGraph::GraphID, double> cluster_sums_of_distances_;

	std::map<GEDGraph::GraphID, double> cluster_radii_;

	double sum_of_distances_;

	Seconds runtime_;

	void set_default_options_();

	void update_medians_();

	void compute_medoid_(const std::vector<GEDGraph::GraphID> & graph_ids, GEDGraph::GraphID medoid_id);

	void update_clusters_();

	double compute_bag_distance_(const std::vector<GEDGraph::GraphID> & cluster_1, const std::vector<GEDGraph::GraphID> & cluster_2) const;

};

}

#include "graph_clustering_heuristic.ipp"

#endif /* MEDIAN_SRC_GRAPH_CLUSTERING_HEURISTICS_HPP_ */
