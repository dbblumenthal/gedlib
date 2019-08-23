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
 * @file compact_mip.hpp
 * @brief ged::CompactMIP class declaration.
 */

#ifndef SRC_METHODS_COMPACT_MIP_HPP_
#define SRC_METHODS_COMPACT_MIP_HPP_


namespace ged {

/*!
 * @brief Mixed integer linear programming formulation of the graph edit distance.
 * @details Implements the compact MIP formulation suggested in:
 * - D. B. Blumenthal and J. Gamper:
 *   &ldquo;On the exact computation of the graph edit distance&rdquo;,
 *   https://doi.org/10.1016/j.patrec.2018.05.002
 *
 * Does not support any options except for the ones supported by ged::MIPBasedMethod.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class CompactMIP : public MIPBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~CompactMIP();

	CompactMIP(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	std::map<NodeMap::Assignment, GRBVar> x_;

	std::map<NodeMap::Assignment, GRBVar> z_;

	// Virtual member functions inherited from MIPBasedMethod.

	virtual void mip_populate_model_(const GEDGraph & g, const GEDGraph & h, GRBModel & model) final;

	virtual void mip_model_to_node_map_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, NodeMap & node_map) final;

	virtual bool mip_model_to_lsape_projection_problem_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, DMatrix & lsape_instance) final;

	// Private helper function.

	char variable_type_() const;

};

}


#endif /* SRC_METHODS_COMPACT_MIP_HPP_ */
