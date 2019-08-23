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
 * @file  branch_tight.hpp
 * @brief ged::BranchTight class declaration.
 */

#ifndef SRC_METHDOS_BRANCH_TIGHT_HPP_
#define SRC_METHDOS_BRANCH_TIGHT_HPP_

namespace ged {

/*!
 * @brief Computes lower and upper bounds for general edit costs.
 * @details Implements the method %BranchTight suggested in:
 * - D. B. Blumenthal and J. Gamper:
 *   &ldquo;Improved lower bounds for graph edit distance&rdquo;,
 *   https:://doi.org/10.1109/TKDE.2017.2772243
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--iterations @<convertible to int greater equal 0@></tt> | maximal number of iterations | @p 20 | https:://doi.org/10.1109/TKDE.2017.2772243 <br> if @p 0, no iteration based termination criterion is used|
 * | <tt>\--time-limit @<convertible to double@></tt> | time limit in seconds | @p 0 | https:://doi.org/10.1109/TKDE.2017.2772243 <br> if less or equal @p 0, no time limit is enforced |
 * | <tt>\--range @<convertible to double@></tt> | range | @p 0 | https:://doi.org/10.1109/TKDE.2017.2772243 <br> if less or equal @p 0, no range based termination criterion is used |
 * | <tt>\--epsilon @<convertible to double@></tt> | range | @p 0 | https:://doi.org/10.1109/TKDE.2017.2772243 <br> if less or equal @p 0, no convergence based termination criterion is used |
 * | <tt>\--regularize NAIVE\|K-FACTOR</tt> | regularization method | @p NAIVE | https:://doi.org/10.1109/TKDE.2017.2772243 |
 * | <tt>\--threads @<convertible to int greater 0@></tt> | number of threads | @p 1 | n.a. |
 * | <tt>\--upper-bound NO\|FIRST\|LAST\|BEST</tt> | upper bound returned by the method | @p BEST | https:://doi.org/10.1109/TKDE.2017.2772243 |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class BranchTight : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:
	virtual ~BranchTight();

	BranchTight(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	enum UpperBoundOption_ {NO, FIRST, LAST, BEST};

	typedef boost::associative_property_map<GEDGraph::NodeNodeMap> MateMap_;

	typedef std::map<std::size_t, std::size_t> SizeTSizeTMap_;

	class Weights_ {

	public:
		Weights_(std::size_t size_master);

		double get_weight(std::size_t row_in_master, std::size_t col_in_master, std::size_t row_sub_in_master, std::size_t col_sub_in_master) const;

		void set_weight(std::size_t row_in_master, std::size_t col_in_master, std::size_t row_sub_in_master, std::size_t col_sub_in_master, double weight);

	private:
		std::size_t size_master_;

		std::vector<double> weights_;

	};

	class SubproblemSolvers_ {

	public:
		SubproblemSolvers_(std::size_t size_master, std::size_t degree);

		LSAPSolver & solver(std::size_t row_in_master, std::size_t col_in_master);

		const LSAPSolver & solver(std::size_t row_in_master, std::size_t col_in_master) const;

		DMatrix & subproblem(std::size_t row_in_master, std::size_t col_in_master);

		const DMatrix & subproblem(std::size_t row_in_master, std::size_t col_in_master) const;

		SizeTSizeTMap_ & rows_subproblem_to_master(std::size_t row_in_master, std::size_t col_in_master);

		std::size_t row_in_master(std::size_t row_in_master, std::size_t col_in_master, std::size_t row_in_subproblem) const;

		SizeTSizeTMap_ & cols_subproblem_to_master(std::size_t row_in_master, std::size_t col_in_master);

		std::size_t col_in_master(std::size_t row_in_master, std::size_t col_in_master, std::size_t col_in_subproblem) const;

		std::size_t get_size() const;

		std::size_t get_degree() const;

		void solve(std::size_t num_threads);

		void clear_solutions();

	private:
		std::size_t size_master_;

		std::size_t degree_;

		std::vector<DMatrix> subproblems_;

		std::vector<LSAPSolver> subproblem_solvers_;

		std::vector<SizeTSizeTMap_> rows_sub_to_master_;

		std::vector<SizeTSizeTMap_> cols_sub_to_master_;
	};

	class KFactor_ {

	public:
		KFactor_(std::size_t num_nodes);

		bool contains_edge(std::size_t id_i, std::size_t id_k) const;

		void add_edge(std::size_t id_i, std::size_t id_k);

		std::size_t num_nodes() const;

		void clear_edges();

	private:
		std::size_t num_nodes_;

		std::vector<bool> is_edge_;
	};

	std::size_t max_itrs_;

	double range_;

	double epsilon_;

	UpperBoundOption_ upper_bound_option_;

	bool naive_regularization_;

	std::size_t num_threads_;

	double time_limit_in_sec_;

	// Member functions inherited from GEDMethod.

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual void ged_set_default_options_() final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ged_valid_options_string_() const;

	// Private helper functions.

	bool termination_criterion_met_(const std::size_t & current_itr, const double & last_improvement, Result & result);

	void init_subproblems_(const GEDGraph & g, const GEDGraph & h, SubproblemSolvers_ & subproblem_solver) const;

	void init_node_costs_(const GEDGraph & g, const GEDGraph & h, DMatrix & node_costs) const;

	void update_master_problem_costs_(const SubproblemSolvers_ & subproblems_solver, const DMatrix & node_costs, DMatrix & master_problem) const;

	void update_subproblem_costs_(const Weights_ & weights, std::size_t degree, SubproblemSolvers_ & subproblems_solver) const;

	void update_weights_(const LSAPSolver & master_problem_solver, std::size_t degree, const SubproblemSolvers_ & subproblems_solver, Weights_ & weights) const;

	std::size_t regularize_(GEDGraph & g, GEDGraph & h) const;

	std::size_t naive_regularize_(GEDGraph & g, GEDGraph & h) const;

	void fill_up_smaller_graph_(GEDGraph & g, GEDGraph & h) const;

	void fill_up_both_graphs_(GEDGraph & g, GEDGraph & h) const;

	void construct_complement_graph_(const GEDGraph & graph, GEDGraph & complement_graph, GEDGraph::NodeNodeMap & complement_to_graph) const;

	void regularize_from_k_factor_(const KFactor_ & k_factor, const GEDGraph::NodeSizeTMap & graph_to_k_factor, GEDGraph & graph) const;

	bool compute_k_factor_(const GEDGraph & complement_graph, LabelID k, const GEDGraph::NodeSizeTMap & complement_graph_to_k_factor, KFactor_ & k_factor) const;

	bool construct_transformed_complement_graph_(const GEDGraph & complement_graph, std::size_t k, GEDGraph & transformed_complement_graph, GEDGraph::NodeNodeMap & transformed_to_original_nodes) const;

};

}

#endif /* SRC_METHDOS_BRANCH_TIGHT_HPP_ */
