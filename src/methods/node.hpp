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
 * @file  node.hpp
 * @brief ged::Node class declaration.
 */

#ifndef SRC_METHODS_NODE_HPP_
#define SRC_METHODS_NODE_HPP_

namespace ged {

/*!
 * \brief Computes lower and upper bounds for general edit costs.
 * \details Implements the method %Node suggested in:
 * - D. Justice and A. Hero:
 *   &ldquo;A binary linear programming formulation of the graph edit distance&rdquo;,
 *   https:://doi.org/10.1109/TPAMI.2006.152
 *
 * Does not support any options except for the ones supported by ged::LSAPEBasedMethod.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Node : public LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Node();

	Node(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	// Member functions inherited from LSAPEBasedMethod.

	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) final;

};

}

#endif /* SRC_METHODS_NODE_HPP_ */
