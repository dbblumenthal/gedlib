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
 * @file  anchor_aware_ged.hpp
 * @brief ged::AnchorAwareGED class declaration.
 */

#ifndef SRC_METHODS_ANCHOR_AWARE_GED_HPP_
#define SRC_METHODS_ANCHOR_AWARE_GED_HPP_

namespace ged {

/*!
 * @brief Computes the exact graph edit distance for general edit costs.
 * @details Implements the following exact algorithm suggested and extends it to non-uniform edit costs:
 * - L. Chang, X. Feng, X. Lin, L. Qin, and W. Zhang:
 *   &ldquo;Efficient graph edit distance computation and verification via anchor-aware lower bound estimation&rdquo;,
 *   https://arxiv.org/abs/1709.06810
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--threads @<convertible to int greater 0@></tt> | number of threads | @p 1 | n.a. |
 * | <tt>\--lsape-model ECBP\|EBP\|FLWC\|FLCC\|FBP\|SFBP\|FBP0</tt> | model for optimally solving LSAPE | @p ECBP | ged::LSAPESolver::Model |
 * | <tt>\--time-limit @<convertible to double@></tt> | time limit in seconds | @p 0 | if less or equal 0, no time limit is enforced |
 * | <tt>\--search-method DFS\|ASTAR</tt> | method for traversing the search tree | @p DFS | https://arxiv.org/abs/1709.06810 <br> if @p DFS, the method computes an upper bound even if the time limit is reached |
 * | <tt>\--lower-bound-method BRANCH\|BRANCH_FAST</tt> | method for computing lower bound for incomplete node map | @p BRANCH_FAST | https://arxiv.org/abs/1709.06810 |
 * | <tt>\--map-root-to-root TRUE\|FALSE</tt> | decide if the roots of the input graphs are mapped to each other | @p FALSE | if @p TRUE, the nodes with ID 0 of the input graphs are mapped to each other |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class AnchorAwareGED : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~AnchorAwareGED();

	AnchorAwareGED(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	enum SearchMethod_ {DFS, ASTAR};

	enum LowerBoundMethod_ {BRANCH, BRANCH_FAST};

	struct Edge_ {
		Edge_(LabelID label, GEDGraph::EdgeID edge_id);

		LabelID label;

		GEDGraph::EdgeID edge_id;

		bool operator<(const Edge_ & rhs) const;
	};

	class SortedEdges_ {
	public:
		SortedEdges_();

		SortedEdges_(const GEDGraph & g);

		void operator=(const SortedEdges_ & rhs);

		const std::vector<Edge_> & get_incident_edges(GEDGraph::NodeID) const;
	private:
		std::map<GEDGraph::NodeID, std::vector<typename AnchorAwareGED<UserNodeLabel, UserEdgeLabel>::Edge_>> sorted_edges_;
	};

	class TreeNode_ {

	public:

		TreeNode_(const GEDGraph & g, const GEDGraph & h, const AnchorAwareGED * exact);

		TreeNode_(const TreeNode_ & tree_node);

		double lower_bound() const;

		void operator=(const TreeNode_ & rhs);

		bool operator<(const TreeNode_ & rhs) const;

		void prepare_for_sibling_generation();

		void prepare_for_child_generation();

		void append_extension(const GEDGraph & g, const GEDGraph & h, const NodeMap & extension);

		void append_next_assignment(const NodeMap & extension);

		void populate_lsape_instance(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance);

		void set_lower_bound_to_leaf(double lower_bound_to_leaf);

		void update_induced_cost(const GEDGraph & g, const GEDGraph & h);

		double induced_cost() const;

		const NodeMap & node_map() const;

		bool is_leaf_node() const;

		void extend_leaf_node(const GEDGraph & g, const GEDGraph & h);

		GEDGraph::NodeID next_unmatched_node_in_g() const;

		GEDGraph::NodeID last_matched_node_in_g() const;

		bool has_unexplored_sibling();

		std::size_t num_unmatched_nodes_in_g() const;

		std::size_t num_unmatched_nodes_in_h() const;

		void update_original_id_of_unmatched_nodes_in_h();

	private:

		const AnchorAwareGED * exact_;

		NodeMap node_map_;

		std::vector<bool> is_matched_node_in_g_;

		std::vector<bool> is_matched_node_in_h_;

		std::vector<bool> is_candidate_in_h_;

		bool dummy_node_is_candidate_in_h_;

		std::vector<std::size_t> original_id_of_unmatched_nodes_in_h_;

		double induced_cost_;

		double lower_bound_to_leaf_;

		std::size_t num_matched_nodes_in_g_;

		std::size_t num_matched_nodes_in_h_;

		double compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const;

		double compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const;

		double compute_branch_fast_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const;

		double compute_branch_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const;
	};

	LSAPESolver::Model lsape_model_;

	SearchMethod_ search_method_;

	LowerBoundMethod_ lower_bound_method_;

	std::size_t num_threads_;

	double time_limit_in_sec_;

	bool map_root_to_root_;

	std::map<GEDGraph::GraphID, SortedEdges_> sorted_edges_;

	NodeMap best_feasible_;

	std::priority_queue<TreeNode_> open_;

	double omega_;

	// Member functions inherited from GEDMethod.

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual void ged_set_default_options_() final;

	virtual std::string ged_valid_options_string_() const final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual void ged_init_() final;

	// Private helper functions.

	void init_graph_(const GEDGraph & graph);

	void generate_best_child_(const GEDGraph & g, const GEDGraph & h, const TreeNode_ & tree_node);

	void generate_best_sibling_(const GEDGraph & g, const GEDGraph & h, const TreeNode_ & tree_node);

	void generate_next_tree_node_(const GEDGraph & g, const GEDGraph & h, TreeNode_ & next_map, bool update_induced_cost, bool update_upper_bound);
};

}

#endif /* SRC_METHODS_ANCHOR_AWARE_GED_HPP_ */
