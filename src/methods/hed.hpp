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
 * @file hed.hpp
 * @brief ged::HED class declaration.
 */

#ifndef SRC_METHODS_HED_HPP_
#define SRC_METHODS_HED_HPP_

namespace ged {

/*!
 * @brief Computes a lower bound for general edit costs.
 * @details Implements the Hausdorff Edit Distance (HED) suggested in:
 * - A. Fischer, C. Y. Suen, V. Frinken, K. Riesen, and H. Bunke:
 *   &ldquo;Approximation of graph edit distance based on Hausdorff matching&rdquo;,
 *   https:://doi.org/10.1016/j.patcog.2014.07.015
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--threads @<convertible to int greater 0@></tt> | number of threads | 1 | can be used by derived classes |
 * | <tt>\--lsape-model ECBP\|EBP\|FLWC\|FLCC\|FBP\|SFBP\|FBP0</tt> | model for optimally solving LSAPE | @p ECBP | ged::LSAPESolver::Model, has no effect unless @p edge-set-distances is set to @p OPTIMAL |
 * | <tt>\--edge-set-distances OPTIMAL\|HED</tt> | method used for computing edge set distances | @p HED | n.a. |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class HED : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~HED();

	HED(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	LSAPESolver::Model lsape_model_;

	std::size_t num_threads_;

	bool use_hed_for_edge_set_distances_;

	// Inherited member functions from GEDMethod.

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ged_valid_options_string_() const final;

	virtual void ged_set_default_options_() final;

	// Helper member functions.

	void populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance) const;

	double compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const;

	double compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const;

	double compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const;
};

}





#endif /* SRC_METHODS_HED_HPP_ */
