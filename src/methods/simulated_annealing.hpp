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
 * @file simulated_annealing.hpp
 * @brief ged::SimulatedAnnealing class declaration.
 */

#ifndef SRC_METHODS_SIMULATED_ANNEALING_HPP_
#define SRC_METHODS_SIMULATED_ANNEALING_HPP_

namespace ged {

/*!
 * @brief Uses LSAPE instances to approximate GED via simulated annealing.
 * @details Implements the simulated annealing method suggested in:
 * - K. Riesen, A. Fischer, H. Bunke:
 *   &ldquo;Improved graph edit distance approximation with simulated annealing&rdquo;,
 *   https:://doi.org/10.1007/978-3-319-58961-9_20,
 *
 * and its extension to general instantiations of LSAPE-GED suggested in:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun:
 *   &ldquo;Designing heuristics for graph edit distance.&rdquo;,
 *   To be submitted to PR.
 *
 *   Supports the following options in addition to the ones supported by ged::MLBasedMethod.
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--threads @<convertible to int greater 0@></tt> | number of threads | 1 | used for populating the LSAPE instance and computing the lower bound |
 * | <tt>\--iterations @<convertible to int greater 0@></tt> | number of iterations | 1000 | https:://doi.org/10.1007/978-3-319-58961-9_20 |
 * | <tt>\--start-probability @<convertible to double between 0 and 1@></tt> | probability of accepting deterioration at start of simulated annealing | 0.8 | https:://doi.org/10.1007/978-3-319-58961-9_20 |
 * | <tt>\--end-probability @<convertible to double between 0 and 1@></tt> | probability of accepting deterioration at end of simulated annealing | 0.01 | https:://doi.org/10.1007/978-3-319-58961-9_20 |
 * | <tt>\--lower-bound-method BRANCH\|BRANCH_FAST\|BRANCH_TIGHT\|NONE</tt> | method for computing lower bound used as a termination criterion | @p NONE | if @p NONE, the lower bound computed by the method for populating the LSAPE instance is used |
 * | <tt>\--lsape-method BIPARTITE\|BRANCH_FAST\|BRANCH_UNIFORM\|BRANCH\|NODE\|RING\|SUBGRAPH\|WALKS\|BIPARTITE_ML\|RING_ML</tt> | method for populating the LSAPE instance | @p BIPARTITE | if @p BIPARTITE, the feature vectors are identical to the ones suggested in https://doi.org/10.1016/j.patrec.2015.10.007 |
 * | <tt>\--lsape-options '[--@<option@> @<arg@>] [...]'</tt> | options string passed to the method used for populating the LSAPE instance | @p '' | ged::Bipartite, ged::BranchFast, ged::BranchUniform, ged::Branch, ged::Node, ged::Ring, ged::Subgraph, ged::Walks |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class SimulatedAnnealing : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~SimulatedAnnealing();

	SimulatedAnnealing(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> * lsape_method_;

	std::string lsape_method_name_;

	std::string lsape_method_options_;

	GEDMethod<UserNodeLabel, UserEdgeLabel> * lower_bound_method_;

	std::string lower_bound_method_name_;

	std::string lower_bound_method_options_;

	std::size_t num_threads_;

	std::size_t num_iterations_;

	double start_probability_;

	double end_probability_;

	// Member functions inherited from GEDMethod.

	virtual void ged_init_() final;

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ged_valid_options_string_() const final;

	virtual void ged_set_default_options_() final;

	// Private helper member functions.

	void generate_candidate_(const GEDGraph & g, const GEDGraph & h, const DMatrix & lsape_instance, const std::vector<std::size_t> & current_order,
			std::vector<std::size_t> & candidate_order, NodeMap & candidate_node_map) const;

};

}

#endif /* SRC_METHODS_SIMULATED_ANNEALING_HPP_ */
