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
 * @file ibd.hpp
 * @brief ged::IBD class declaration.
 */

#ifndef SRC_EDIT_COSTS_IBD_HPP_
#define SRC_EDIT_COSTS_IBD_HPP_

#include "edit_costs.hpp"
#include "../util/misc.hpp"

namespace ged {

/*!
 * @brief Edit cost functions for IBD graphs.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class IBD : public EditCosts<UserNodeLabel, UserEdgeLabel> {
public:

	virtual ~IBD();

	/*!
	 * @brief Constructor.
	 * @param[in] otu_distances Path to CSV file with the distances between the OTUs.
	 * @param[in] alpha Controls importance of node edit operations vs. importance of edge edit operations.
	 */
	IBD(const std::string & otu_distances, double alpha = 0.5);

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

	double alpha_;

	std::vector<std::size_t> otus_;

};

}

#include "ibd.ipp"



#endif /* SRC_EDIT_COSTS_IBD_HPP_ */
