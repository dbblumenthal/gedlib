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
 * @file  hybrid.hpp
 * @brief ged::Hybrid class declaration.
 */

#ifndef SRC_METHODS_HYBRID_HPP_
#define SRC_METHODS_HYBRID_HPP_

namespace ged {

/*!
 * @brief Computes a lower bound for uniform edit costs.
 * @details Implements the method %Hybrid suggested in:
 * - W. Zheng, L. Zou, X. Lian, D. Wang, and D. Zhao:
 *   &ldquo;Efficient graph similarity search over large graph databases&rdquo;,
 *   https://doi.org/10.1109/TKDE.2014.2349924
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--threads @<convertible to int greater 0@></tt> | number of threads | @p 1 | n.a. |
 * | <tt>\--lsape-model ECBP\|EBP\|FLWC\|FLCC\|FBP\|SFBP\|FBP0</tt> | model for optimally solving LSAPE | @p ECBP | ged::LSAPESolver::Model |
 * | <tt>\--time-limit @<convertible to double@></tt> | time limit in seconds | @p 0 | if less or equal @p 0, no time limit is enforced |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Hybrid : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Hybrid();

	Hybrid(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	typedef typename std::vector<typename Partition<UserNodeLabel, UserEdgeLabel>::Substruct_>::const_iterator SubstructItr_;

	LSAPESolver::Model lsape_model_;

	std::size_t num_threads_;

	double time_limit_in_sec_;

	// Member functions inherited from GEDMethod.

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual void ged_set_default_options_() final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ged_valid_options_string_() const final;

	// Private helper member functions.
	bool branch_uniform_dfs_(const Timer & timer, const std::string & options, GEDGraph & g, GEDGraph & h, SubstructItr_ current_substruct, SubstructItr_ end_substructs, double & lower_bound);

	std::string to_string_(LSAPESolver::Model lsape_model);

};

}

#endif /* SRC_METHODS_HYBRID_HPP_ */
