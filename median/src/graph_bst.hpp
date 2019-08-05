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

#include "graph_clustering_heuristic.hpp"

namespace ged {

template<class UserNodeID, class UserNodeLabel, class UserEdgeLabel>
class GraphBST {

public:

	GraphBST(GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env, GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel> * clustering_heuristic);

	~GraphBST();

	void init(std::vector<GEDGraph::GraphID> graph_ids, std::size_t max_cluster_size);

	void process_range_query(GEDGraph::GraphID query_graph_id, double threshold);

	void save(const std::string & focal_graph_dir, const std::string & bst_file_name) const;

	void load(const std::string & data_graph_dir, const std::string & focal_graph_dir, const std::string & bst_file_name) const;

	void set_lower_bound_method(Options::GEDMethod lower_bound_method, const std::string & lower_bound_options = "");

	void set_upper_bound_method(Options::GEDMethod upper_bound_method, const std::string & upper_bound_options = "");

	const std::vector<GEDGraph::GraphID> & get_verified_graph_ids() const;

	const std::vector<GEDGraph::GraphID> & get_filtered_graph_ids() const;

	const std::vector<GEDGraph::GraphID> & get_undecided_graph_ids() const;

	double get_init_time() const;

	double get_query_time() const;

private:

	struct BSTNode_ {

		BSTNode_();

		~BSTNode_();

		GEDGraph::GraphID focal_graph_id;

		double radius;

		std::vector<GEDGraph::GraphID> cluster;

		GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::BSTNode_ * child_1;

		GraphBST<UserNodeID, UserNodeLabel, UserEdgeLabel>::BSTNode_ * child_2;
	};

	GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> * ged_env_;

	GraphClusteringHeuristic<UserNodeID, UserNodeLabel, UserEdgeLabel> * clustering_heuristic_;

	Options::GEDMethod lower_bound_method_;

	std::string lower_bound_options_;

	Options::GEDMethod upper_bound_method_;

	std::string upper_bound_options_;

	BSTNode_ * root_;

	std::vector<GEDGraph::GraphID> focal_graph_ids_;

	std::vector<GEDGraph::GraphID> verified_graph_ids_;

	std::vector<GEDGraph::GraphID> filtered_graph_ids_;

	std::vector<GEDGraph::GraphID> undecided_graph_ids_;

	Seconds init_time_;

	Seconds query_time_;

	void init_bst_node_(BSTNode_ * bst_node, GEDGraph::GraphID focal_graph_id, double radius, const std::vector<GEDGraph::GraphID> & cluster);

	std::size_t split_(std::size_t max_cluster_size, const BSTNode_ * node, std::size_t pos_next_focal_graph_id);

	void check_(GEDGraph::GraphID query_graph_id, double threshold, const BSTNode_ * node);

	void serialize_(const BSTNode_ * node, const std::string & focal_graph_dir, std::ofstream & file);

	void de_serialize_(const BSTNode_ * node, const std::string & data_graph_dir, const std::string & focal_graph_dir, std::ofstream & file);

};

}

#include "graph_bst.ipp"

#endif /* MEDIAN_SRC_GRAPH_BST_HPP_ */
