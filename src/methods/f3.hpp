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
 * @file f3.hpp
 * @brief ged::F3 class declaration.
 */

#ifndef SRC_METHODS_F3_HPP_
#define SRC_METHODS_F3_HPP_

namespace ged {

/*!
 * @brief Mixed integer linear programming formulation of the graph edit distance.
 * @details Implements the MIP formulation suggested in:
 * - M. Darwiche, R. Raveaux, D. Conte, V. T'Kindt:
 *   &ldquo;Graph edit distance in the exact context&rdquo;,
 *  https://doi.org/10.1007/978-3-319-97785-0_29
 *
 * Does not support any options except for the ones supported by ged::MIPBasedMethod.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class F3 : public MIPBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~F3();

	F3(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	std::map<std::pair<GEDGraph::NodeID, GEDGraph::NodeID>, GRBVar> x_;

	std::map<std::pair<GEDGraph::EdgeID, std::pair<GEDGraph::NodeID, GEDGraph::NodeID>>, GRBVar> y_;

	double constant_;

	// Virtual member functions inherited from MIPBasedMethod.

	virtual void mip_populate_model_(const GEDGraph & g, const GEDGraph & h, GRBModel & model) final;

	virtual double mip_objective_constant_(const GEDGraph & g, const GEDGraph & h) final;

	virtual void mip_model_to_node_map_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, NodeMap & node_map) final;

	// Private helper function.

	char variable_type_() const;

};

}


#endif /* SRC_METHODS_F3_HPP_ */
