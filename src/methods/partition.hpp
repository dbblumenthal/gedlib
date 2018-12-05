/***************************************************************************
*                                                                          *
*   Copyright (C) 2018 by David B. Blumenthal                              *
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
 *  @file  partition.hpp
 *  @brief Partition class declaration.
 */

#ifndef SRC_METHODS_PARTITION_HPP_
#define SRC_METHODS_PARTITION_HPP_

namespace ged {

/*!
 * @brief Computes a lower bound for uniform edit costs.
 * @details Implements the method %Partition suggested in:
 * - W. Zheng, L. Zou, X. Lian, D. Wang, and D. Zhao:
 *   &ldquo;Efficient graph similarity search over large graph databases&rdquo;,
 *   https:://doi.org/10.1109/TKDE.2014.2349924
 *
 * This method does not support any options.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Partition : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

	template<typename HybridUserNodeLabel, typename HybridUserEdgeLabel>
	friend class Hybrid;

public:

	virtual ~Partition();

	Partition(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	enum SubstructType_ {NODE, EDGE, NODE_EDGE, NODE_EDGE_NODE, NODE_EDGE_EDGE};

	struct Substruct_ {
		Substruct_(SubstructType_ type, LabelID label_1, LabelID label_2 = dummy_label(), LabelID label_3 = dummy_label(), GEDGraph::NodeID node_1 = GEDGraph::dummy_node(),
				GEDGraph::EdgeID edge_1 = GEDGraph::dummy_edge(), GEDGraph::NodeID node_2 = GEDGraph::dummy_node(), GEDGraph::EdgeID edge_2 = GEDGraph::dummy_edge());

		SubstructType_ type;

		LabelID node_label_1;

		LabelID node_label_2;

		LabelID edge_label_1;

		LabelID edge_label_2;

		GEDGraph::NodeID node_1;

		GEDGraph::NodeID node_2;

		GEDGraph::EdgeID edge_1;

		GEDGraph::EdgeID edge_2;

		bool operator<(const Substruct_ & rhs) const;
	};

	typedef std::map<Substruct_, bool> SubstructMap_;

	std::map<GEDGraph::GraphID, SubstructMap_> substruct_maps_;

	std::vector<Substruct_> unmatched_substructs_;

	// Member functions inherited from GEDMethod.

	virtual void ged_init_() final;

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	// Private helper member functions.

	void init_graph_(const GEDGraph & graph);

	void init_graphs_(const GEDGraph & g, const GEDGraph & h);

	const std::vector<Substruct_> & get_unmatched_substructs_() const;

	void check_node_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::NodeID, bool> & is_deleted_node);

	void check_edge_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::EdgeID, bool> & is_deleted_edge);

	void check_node_edge_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::NodeID, bool> & is_deleted_node, std::map<GEDGraph::EdgeID, bool> & is_deleted_edge);

	void check_node_edge_node_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::NodeID, bool> & is_deleted_node, std::map<GEDGraph::EdgeID, bool> & is_deleted_edge);

	void check_node_edge_edge_subtructs_(const GEDGraph & h, const SubstructMap_ & is_substruct_in_g, std::map<GEDGraph::NodeID, bool> & is_deleted_node, std::map<GEDGraph::EdgeID, bool> & is_deleted_edge);

	bool contains_substruct_(const GEDGraph & g, const Substruct_ & substruct) const;

	bool contains_node_substruct_(const GEDGraph & g, const Substruct_ & substruct) const;

	bool contains_edge_substruct_(const GEDGraph & g, const Substruct_ & substruct) const;

	bool contains_node_edge_substruct_(const GEDGraph & g, const Substruct_ & substruct) const;

	bool contains_node_edge_node_substruct_(const GEDGraph & g, const Substruct_ & substruct) const;

	bool contains_node_edge_edge_substruct_(const GEDGraph & g, const Substruct_ & substruct) const;

};

}

#endif /* SRC_METHODS_PARTITION_HPP_ */
