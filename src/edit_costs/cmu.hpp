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
 * @file cmu.hpp
 * @brief ged::CMU class declaration.
 */

#ifndef SRC_EDIT_COSTS_CMU_HPP_
#define SRC_EDIT_COSTS_CMU_HPP_

#include "edit_costs.hpp"

namespace ged {

/*!
 * @brief Edit costs for graphs contain in CMU dataset.
 * @details The graphs contained in the CMU dataset represent images of houses from different viewpoints.
 * Each graph has 30 nodes that are attributed with Euclidian coordinates (named "x" and "y").
 * An edge e is attributed with the integer-valued Euclidean distance d(e) between the attributes of its endpoints (named "dist").
 * The CMU dataset is part of the graph data repository for graph edit distance which can be downloaded from http://www.rfai.li.univ-tours.fr/PublicData/GDR4GED/home.html:
 * - Z. Abu-Aisheh, R. Raveaux, J.-Y. Ramel:
 *   &ldquo;A graph database repository and performance evaluation metrics for graph edit distance&rdquo;,
 *   https://doi.org/10.1007/978-3-319-18224-7_14
 *
 * Implements the edit costs suggested in:
 * - Z. Abu-Aisheh, B. Gaüzère, S. Bougleux, J.-Y. Ramel, L. Brun, R. Raveaux, P. Héroux, and S. Adam.
 *   &ldquo;Graph edit distance contest 2016: Results and future challenges&rdquo;,
 *   https://doi.org/10.1016/j.patrec.2017.10.007
 * - Z. Abu-Aisheh, R. Raveaux, J.-Y. Ramel:
 *   &ldquo;A graph database repository and performance evaluation metrics for graph edit distance&rdquo;,
 *   https://doi.org/10.1007/978-3-319-18224-7_14
 */
template<class UserNodeLabel, class UserEdgeLabel>
class CMU : public EditCosts<UserNodeLabel, UserEdgeLabel> {
public:

	virtual ~CMU();

	/*!
	 * @brief Constructor.
	 * @param[in] node_ins_del_cost Cost for deleting or inserting nodes.
	 * @param[in] alpha Importance of node edit operations vs. importance of edge edit operations.
	 * @note Calling the constructor with the default arguments constructs the edit costs suggested in https://doi.org/10.1016/j.patrec.2017.10.007 and https://doi.org/10.1007/978-3-319-18224-7_14.
	 */
	CMU(double node_ins_del_cost = 100000, double alpha = 0.5);

	virtual double node_ins_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_del_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const final;

	virtual void vectorize_node_label(const UserNodeLabel & node_label, std::vector<double> & vector_representation) const final;

	virtual double edge_ins_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_del_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_rel_cost_fun(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const final;

private:

	double node_ins_del_cost_;

	double alpha_;
};

}

#include "cmu.ipp"

#endif /* SRC_EDIT_COSTS_CMU_HPP_ */
