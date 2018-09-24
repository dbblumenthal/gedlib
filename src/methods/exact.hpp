/*!
 * @file  exact.hpp
 * @brief ged::Exact class declaration.
 */

#ifndef SRC_METHODS_EXACT_HPP_
#define SRC_METHODS_EXACT_HPP_

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
class Exact : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Exact();

	Exact(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

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
		std::map<GEDGraph::NodeID, std::vector<typename Exact<UserNodeLabel, UserEdgeLabel>::Edge_>> sorted_edges_;
	};

	struct NodeMap_ {
		NodeMap_(const GEDGraph & g, const GEDGraph & h, const Exact * exact);

		NodeMap_(const Exact * exact);

		NodeMap_(const NodeMap_ & node_map);

		const Exact * exact;

		NodeMap matching;

		std::map<GEDGraph::NodeID, bool> is_matched_node_in_g;

		std::map<GEDGraph::NodeID, bool> is_matched_node_in_h;

		std::map<GEDGraph::NodeID, bool> is_candidate_in_h;

		double induced_cost;

		double lower_bound_to_leaf;

		std::size_t num_matched_nodes_in_g;

		std::size_t num_matched_nodes_in_h;

		double lower_bound() const;

		void operator=(const NodeMap_ & rhs);

		bool operator<(const NodeMap_ & rhs) const;

		GEDGraph::NodeID next_unmatched_node_in_g(const GEDGraph & g) const;

		GEDGraph::NodeID last_matched_node_in_g(const GEDGraph & g) const;

		void reset_is_candidate_in_h();

		bool candidates_left();

		void print();
	};

	LSAPESolver::Model lsape_model_;

	SearchMethod_ search_method_;

	LowerBoundMethod_ lower_bound_method_;

	std::size_t num_threads_;

	double time_limit_in_sec_;

	bool map_root_to_root_;

	std::map<GEDGraph::GraphID, SortedEdges_> sorted_edges_;

	NodeMap_ best_feasible_;

	std::priority_queue<NodeMap_> open_;

	double omega_;

	// Member functions inherited from GEDMethod.

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual void ged_set_default_options_() final;

	virtual std::string ged_valid_options_string_() const final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual void ged_init_() final;

	// Private helper functions.

	void init_graph_(const GEDGraph & graph);

	void extend_half_complete_node_map_(const GEDGraph & g, const GEDGraph & h, NodeMap_ & node_map) const;

	void generate_best_child_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & node_map);

	void generate_best_sibling_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & node_map);

	void generate_next_map_(const GEDGraph & g, const GEDGraph & h, NodeMap_ & next_map, bool update_induced_cost, bool update_upper_bound);

	void init_indices_(const NodeMap_ & node_map, GEDGraph::SizeTNodeMap & g_ids_to_nodes, GEDGraph::SizeTNodeMap & h_ids_to_nodes) const;

	void init_master_problem_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & node_map, const GEDGraph::SizeTNodeMap & g_ids_to_nodes, const GEDGraph::SizeTNodeMap &  h_ids_to_nodes, DMatrix & master_problem) const;

	double compute_insertion_cost_(const GEDGraph & h, const NodeMap_ & node_map, GEDGraph::NodeID k) const;

	double compute_deletion_cost_(const GEDGraph & g, const NodeMap_ & node_map, GEDGraph::NodeID i) const;

	double compute_branch_substitution_cost_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & node_map, GEDGraph::NodeID i, GEDGraph::NodeID k) const;

	double compute_branch_fast_substitution_cost_(const GEDGraph & g, const GEDGraph & h, const NodeMap_ & node_map, GEDGraph::NodeID i, GEDGraph::NodeID k) const;

	void append_extension_(const GEDGraph & g, const GEDGraph & h, const NodeMap & extension, NodeMap_ & node_map);

	void update_induced_cost_(const GEDGraph & g, const GEDGraph & h, NodeMap_ & child_map) const;
};

}

#endif /* SRC_METHODS_EXACT_HPP_ */
