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
 * @file  lsape_based_method.hpp
 * @brief ged::LSAPEBasedMethod class declaration.
 */

#ifndef SRC_METHODS_LSAPE_BASED_METHOD_HPP_
#define SRC_METHODS_LSAPE_BASED_METHOD_HPP_

namespace ged {

/*!
 * @brief Abstract class for methods that use lossy transformations to LSAPE for approximating the graph edit distance.
 * @details Implements the paradigm LSAPE-GED described in:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, and L. Brun:
 *   &ldquo;%Ring based approximation of graph edit distance&rdquo;,
 *   S+SSPR 2018,
 *
 * and the extension of LSAPE-GED which uses node centrality measures:
 * - K. Riesen, H. Bunke, and A. Fischer:
 *   &ldquo;Improving graph edit distance approximation by centrality measures&rdquo;,
 *   https:://doi.org/10.1016/j.patcog.2014.11.002
 *
 * All derived classes support the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--threads @<convertible to int greater 0@></tt> | number of threads | 1 | can be used by derived classes |
 * | <tt>\--lsape-model ECBP\|EBP\|FLWC\|FLCC\|FBP\|SFBP\|FBP0</tt> | model for optimally solving LSAPE | @p ECBP | ged::LSAPESolver::Model |
 * | <tt>\--greedy-method BASIC\|REFINED\|LOSS\|BASIC_SORT\|INT_BASIC_SORT</tt> | method for greedily solving LSAPE  | @p BASIC | ged::LSAPESolver::GreedyMethod |
 * | <tt>\--optimal TRUE\|FALSE</tt> | solve optimally or greedily | @p TRUE | if @p TRUE, the option @p \--greedy-method has no effect <br> if @p FALSE, the option @p \--lsape-method has no effect |
 * | <tt>\--centrality-method NONE\|DEGREE\|EIGENVECTOR\|PAGERANK</tt> | node centrality method | @p NONE | https:://doi.org/10.1016/j.patcog.2014.11.002 <br> if @p NONE, the option @p \--centrality-weight has no effect|
 * | <tt>\--centrality-weight @<convertible to double between 0 and 1@></tt> | weight of node centralities | 0.7 | https:://doi.org/10.1016/j.patcog.2014.11.002 |
 * | <tt>\--max-num-solutions ALL\|@<convertible to int greater 0@></tt> | maximal number of solutions | 1 | the actual number of computed solutions might be smaller than the specified value |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class LSAPEBasedMethod : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	// Constructor and destructor.

	/*!
	 * @brief Pure virtual destructor.
	 * @note Must be implemented by derived classes.
	 */
	virtual ~LSAPEBasedMethod() = 0;

	/*!
	 * @brief Constructor.
	 * @param[in] ged_data The instance on which the method should be run.
	 */
	LSAPEBasedMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

	// Extension of interface provided by GEDMethod.

	/*!
	 * @brief Runs the method with options specified by set_options() and provides access to constructed LSAPE instance.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[out] result Result variable.
	 * @param[out] lsape_instance LSAPE instance.
	 */
	void populate_instance_and_run_as_util(const GEDGraph & g, const GEDGraph & h, Result & result, DMatrix & lsape_instance);

	/*!
	 * @brief Populates the LSAPE instance.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[out] lsape_instance LSAPE instance.
	 */
	void populate_instance(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance);

protected:

	/*!
	 * @brief Specifies model for optimal LSAPE solver.
	 */
	LSAPESolver::Model lsape_model_;

	/*!
	 * @brief Specifies method for greedy LSAPE solver.
	 */
	LSAPESolver::GreedyMethod greedy_method_;

	/*!
	 * @brief Specifies the enumeration method for optimal LSAPE solver.
	 */
	LSAPESolver::EnumerationMethod enumeration_method_;

	/*!
	 * @brief Flag that should be set to @p true if and only if the method computes a lower bound.
	 */
	bool compute_lower_bound_;

	/*!
	 * @brief Flag that equals @p true if an optimal LSAPE solver is used and @p false if a greedy method is employed.
	 */
	bool solve_optimally_;

	/*!
	 * @brief The number of threads to be used.
	 */
	std::size_t num_threads_;

private:

	enum CentralityMethod_ {NONE, DEGREE, EIGENVECTOR, PAGERANK};

	CentralityMethod_ centrality_method_;

	double centrality_weight_;

	std::map<GEDGraph::GraphID, std::vector<double>> centralities_;

	int max_num_solutions_;

	// Member functions inherited from GEDMethod.

	virtual void ged_init_() final;

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ged_valid_options_string_() const final;

	virtual void ged_set_default_options_() final;

	// Private helper member functions.

	void init_graph_(const GEDGraph & graph);

	void init_centralities_(const GEDGraph & graph);

	void add_centralities_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem);

	void compute_eigenvector_with_largest_eigenvalue_(const DMatrix & symmetric_matrix, std::vector<double> & eigenvector, double & eigenvalue);

	// Virtual member functions to be overridden by derived classes.

	/*!
	 * @brief Initializes the method after initializing the global variables for the graphs.
	 * @note Must be overridden by derived classes of ged::LSAPEBasedMethod that require custom initialization.
	 */
	virtual void lsape_init_();

	/*!
	 * @brief Parses one option that is not among the ones shared by all derived classes of ged::LSAPEBasedMethod.
	 * @param[in] option The name of the option.
	 * @param[in] arg The argument of the option.
	 * @return Returns true if @p option is a valid option name for the method and false otherwise.
	 * @note Must be overridden by derived classes of ged::LSAPEBasedMethod that have options that are not among the ones shared by all derived classes of ged::LSAPEBasedMethod.
	 */
	virtual bool lsape_parse_option_(const std::string & option, const std::string & arg);

	/*!
	 * @brief Returns string of all valid options that are not among the ones shared by all derived classes of ged::LSAPEBasedMethod.
	 * @return String of the form @"[--@<option@> @<arg@>] [...]@".
	 * @note Must be overridden by derived classes of ged::LSAPEBasedMethod that have options that are not among the ones shared by all derived classes of ged::LSAPEBasedMethod.
	 */
	virtual std::string lsape_valid_options_string_() const;

	/*!
	 * @brief Sets all options that are not among the ones shared by all derived classes of ged::LSAPEBasedMethod to default values.
	 * @note Must be overridden by derived classes of ged::LSAPEBasedMethod that have options that are not among the ones shared by all derived classes of ged::LSAPEBasedMethod.
	 */
	virtual void lsape_set_default_options_();

	/*!
	 * @brief Populates the LSAPE instance.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[out] lsape_instance LSAPE instance of size (n + 1) x (m + 1), where n and m are the number of nodes in @p g and @p h. The last row and the last column represent insertion and deletion.
	 * @note Must be overridden by derived classes of ged::LSAPEBasedMethod.
	 */
	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance);

	/*!
	 * @brief Initializes global variables for one graph.
	 * @param[in] graph Graph for which the global variables have to be initialized.
	 * @note Must be overridden by derived classes of ged::LSAPEBasedMethod that require to initialize custom global variables.
	 */
	virtual void lsape_init_graph_(const GEDGraph & graph);

	/*!
	 * @brief Initializes the method at runtime or during initialization before initializing the global variables for the graphs.
	 * @param[in] called_at_runtime Equals @p true if called at runtime and @p false if called during initialization.
	 * @brief Must be overridden by derived classes of ged::LSAPEBasedMethod that require default initialization at runtime before initializing the global variables for the graphs.
	 */
	virtual void lsape_pre_graph_init_(bool called_at_runtime);

	/*!
	 * @brief Default initializes the method at runtime after initializing the global variables for the graphs.
	 * @note Must be overridden by derived classes of ged::LSAPEBasedMethod that require default initialization at runtime after initializing the global variables for the graphs.
	 */
	virtual void lsape_default_post_graph_init_();

	/*!
	 * @brief Returns scaling factor for lower bound.
	 * @param[in] g g Input graph.
	 * @param[in] h h Input graph.
	 * @return The factor by which the optimal LSAPE solution has to be scaled in order to arrive at a valid lower bound.
	 * @note Must be overridden by derived classes of ged::LSAPEBasedMethod that require scaling of the optimal LSAPE solution for computing their lower bounds.
	 */
	virtual double lsape_lower_bound_scaling_factor_(const GEDGraph & g, const GEDGraph & h);

};

}

#endif /* SRC_METHODS_LSAPE_BASED_METHOD_HPP_ */
