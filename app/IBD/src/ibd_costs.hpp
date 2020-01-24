/***************************************************************************
*                                                                          *
*   Copyright (C) 2020 by David B. Blumenthal                              *
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
 * @file ibd_costs.hpp
 * @brief ged::IBDCosts class declaration.
 */

#ifndef APP_IBD_SRC_IBD_COSTS_HPP_
#define APP_IBD_SRC_IBD_COSTS_HPP_

#include "../../../src/edit_costs/edit_costs.hpp"
#include "../../../src/util/misc.hpp"

using namespace ged;

/*!
 * @brief Edit cost functions for IBD graphs.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class IBDCosts : public EditCosts<UserNodeLabel, UserEdgeLabel> {
public:

	virtual ~IBDCosts();

	/*!
	 * @brief Constructor.
	 * @param[in] distance_matrix Path to distance matrix for node relabeling costs.
	 * @param[in] ins_del_factor Importance of insertions and deletions w.r.t. importance of substitutions.
	 * @param[in] alpha Importance of node edit operations vs. importance of edge edit operations.
	 */
	IBDCosts(const std::string & distance_matrix, double ins_del_factor = 1, double alpha = 0.5);

	virtual double node_ins_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_del_cost_fun(const UserNodeLabel & node_label) const final;

	virtual double node_rel_cost_fun(const UserNodeLabel & node_label_1, const UserNodeLabel & node_label_2) const final;

	virtual UserNodeLabel median_node_label(const std::vector<UserNodeLabel> & node_labels) const final;

	virtual double edge_ins_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_del_cost_fun(const UserEdgeLabel & edge_label) const final;

	virtual double edge_rel_cost_fun(const UserEdgeLabel & edge_label_1, const UserEdgeLabel & edge_label_2) const final;

	virtual UserEdgeLabel median_edge_label(const std::vector<UserEdgeLabel> & edge_labels) const final;

private:

	DMatrix node_rel_costs_;

	double max_node_rel_cost_;

	std::vector<std::size_t> otu_to_index_;

	std::vector<std::size_t> index_to_otu_;

	double max_edge_rel_cost_;

	double alpha_;

	double ins_del_factor_;
};


#include "ibd_costs.ipp"



#endif /* APP_IBD_SRC_IBD_COSTS_HPP_ */
