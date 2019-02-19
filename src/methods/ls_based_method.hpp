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
 * @file ls_based_method.hpp
 * @brief ged::LSBasedMethod class declaration.
 */

#ifndef SRC_METHODS_LS_BASED_METHOD_HPP_
#define SRC_METHODS_LS_BASED_METHOD_HPP_

namespace ged {

/*!
 * @brief Abstract class for methods that use variants of local search for upper bounding the graph edit distance.
 * @details All derived methods support the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--initialization-method BIPARTITE_ML\|BIPARTITE\|BRANCH_FAST\|BRANCH_UNIFORM\|BRANCH\|NODE\|RING_ML\|RING\|SUBGRAPH\|WALKS\|RANDOM</tt> | method for computing the initial solution | @p RANDOM | if @p RANDOM, the option @p \--initialization-options has no effect |
 * | <tt>\--initialization-options '[--@<option@> @<arg@>] [...]'</tt> | options string passed to the initialization method | @p '' | ged::BipartiteML, ged::Bipartite, ged::BranchFast, ged::BranchUniform, ged::Branch, ged::Node ged::RingML, ged::Ring, ged::Subgraph, ged::Walks |
 * | <tt>\--lower-bound-method BRANCH\|BRANCH_FAST\|BRANCH_TIGHT\|NONE</tt> | method for computing lower bound used as a termination criterion | @p NONE | if @p NONE, the lower bound computed by the initialization method is used |
 * | <tt>\--random-substitution-ratio @<convertible to double between 0 and 1@></tt> | ratio of node substitutions in randomly constructed initial solutions | @p 1 | n.a |
 * | <tt>\--initial-solutions @<convertible to int greater 0@></tt> | number of initial solutions | @p 1 | if greater @p 1 and if a non-random initialization method is chosen, those initial solution that cannot be computed non-randomly because the initialization method does not yield enough solutions are computed randomly |
 * | <tt>\--ratio-runs-from-initial-solutions @<convertible to double greater 0 and smaller equal 1@></tt> | determines number of runs from initial solutions | @p 1 | if set to @p r and @p \--initial-solutions is set to <tt>k</tt>, all remaining initial solutions are discarded as soon as @p r * @p k runs have been completed |
 * | <tt>\--threads @<convertible to int greater 0@></tt> | number of threads | @p 1 | used for initializing the initialization method and for parallelly running local searches from several initial solutions |
 * | <tt>\--num-randpost-loops @<convertible to int greater equal 0@></tt> | number of loops of RANDPOST algorithm | @p 0 | https://doi.org/10.1007/978-3-319-97785-0_44 |
 * | <tt>\--max-randpost-retrials @<convertible to int greater equal 0@></tt> | number of times the RANDPOST algorithm flattens the probability distribution if it encounters already converged solutions | @p 10 | https://doi.org/10.1007/978-3-319-97785-0_44 |
 * | <tt>\--randpost-penalty @<convertible to double between 0 and 1@></tt> | if set value close to @p 1, expensive solutions count less for constructing the counts matrix | @p 0 | n.a. |
 * | <tt>\--randpost-decay @<convertible to double between 0 and 1@></tt> | if set value close to @p 0, previous iterations contribute less when constructing the counts matrix | @p 1 | n.a. |
 * | <tt>\--log @<file name></tt> | name of log-file for RANDPOST | @p "" | if not empty, a log-file for RANDPOST is created |
 * | <tt>\--randomness REAL\|PSEUDO</tt> | use real randomness or pseudo randomness | @p REAL | n.a. |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class LSBasedMethod : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	/*!
	 * @brief Pure virtual destructor.
	 * @note Must be implemented by derived classes.
	 */
	virtual ~LSBasedMethod() = 0;

	/*!
	 * @brief Constructor.
	 * @param[in] ged_data The instance on which the method should be run.
	 */
	LSBasedMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

protected:

	/*!
	 * @brief The number of threads to be used.
	 */
	std::size_t num_threads_;

private:

	LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> * initialization_method_;

	std::string initialization_options_;

	GEDMethod<UserNodeLabel, UserEdgeLabel> * lower_bound_method_;

	std::string lower_bound_method_options_;

	double random_substitution_ratio_;

	std::size_t num_initial_solutions_;

	double ratio_runs_from_initial_solutions_;

	std::size_t num_randpost_loops_;

	std::size_t max_randpost_retrials_;

	double randpost_penalty_;

	double randpost_decay_;

	std::string logfile_name_;

	bool use_real_randomness_;

	// Member functions inherited from GEDMethod.

	virtual void ged_init_() final;

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ged_valid_options_string_() const final;

	virtual void ged_set_default_options_() final;

	// Private helper member functions.

	std::size_t num_runs_from_initial_solutions_() const;

	void generate_initial_node_maps_(const GEDGraph & g, const GEDGraph & h, std::vector<NodeMap> & initial_node_maps, Result & result);

	void generate_random_initial_node_maps_(const GEDGraph & g, const GEDGraph & h, std::vector<NodeMap> & initial_node_maps);

	void generate_lsape_based_initial_node_maps_(const GEDGraph & g, const GEDGraph & h, std::vector<NodeMap> & initial_node_maps, Result & result);

	double update_counts_matrix_and_visited_node_maps_(const std::vector<NodeMap> & result_node_maps, const std::vector<bool> & is_converged_node_map, const double & upper_bound,
			const double & lower_bound, std::vector<NodeMap> & visited_node_maps, std::size_t loop, std::vector<std::vector<double>> & counts_matrix) const;

	void generate_node_maps_from_counts_matrix_(const GEDGraph & g, const GEDGraph & h,const std::vector<std::vector<double>> & counts_matrix, std::vector<NodeMap> & visited_node_maps, std::vector<NodeMap> & initial_node_maps) const;

	// Virtual member functions to be overridden by derived classes.

	/*!
	 * @brief Runs the local search from an initial node map.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[in] lower_bound A lower bound for GED which can be used as a termination criterion.
	 * @param[in] initial_node_map Initial node map.
	 * @param[out] output_node_map Node map constructed by local search.
	 * @note Must be overridden by derived classes of ged::LSBasedMethod.
	 */
	virtual void ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map);

	/*!
	 * @brief Initializes the method.
	 * @note Must be overridden by derived classes of ged::LSBasedMethod that require initialization.
	 */
	virtual void ls_init_();

	/*!
	 * @brief Initializes the method for a run between two graphs.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @note Must be overridden by derived classes of ged::LSBasedMethod that require initialization at runtime.
	 */
	virtual void ls_runtime_init_(const GEDGraph & g, const GEDGraph & h);

	/*!
	 * @brief Parses one option that is not among the ones shared by all derived classes of ged::LSBasedMethod.
	 * @param[in] option The name of the option.
	 * @param[in] arg The argument of the option.
	 * @return Returns true if @p option is a valid option name for the method and false otherwise.
	 * @note Must be overridden by derived classes of ged::LSBasedMethod that have options that are not among the ones shared by all derived classes of ged::LSBasedMethod.
	 */
	virtual bool ls_parse_option_(const std::string & option, const std::string & arg);

	/*!
	 * @brief Returns string of all valid options that are not among the ones shared by all derived classes of ged::LSBasedMethod.
	 * @return String of the form @"[--@<option@> @<arg@>] [...]@".
	 * @note Must be overridden by derived classes of ged::LSBasedMethod that have options that are not among the ones shared by all derived classes of ged::LSBasedMethod.
	 */
	virtual std::string ls_valid_options_string_() const;

	/*!
	 * @brief Sets all options that are not among the ones shared by all derived classes of ged::LSBasedMethod to default values.
	 * @note Must be overridden by derived classes of ged::LSBasedMethod that have options that are not among the ones shared by all derived classes of ged::LSBasedMethod.
	 */
	virtual void ls_set_default_options_();

};

}

#endif /* SRC_METHODS_LS_BASED_METHOD_HPP_ */
