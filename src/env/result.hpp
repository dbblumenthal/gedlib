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
 * @file result.hpp
 * @brief ged::Result struct declaration.
 */

#ifndef SRC_ENV_RESULT_HPP_
#define SRC_ENV_RESULT_HPP_

#include "common_types.hpp"
#include "node_map.hpp"

namespace ged {

/*!
 * @brief A wrapper structure for the result of calls to ged::GEDMethod::run_as_util() and ged::GEDMethod::ged_run_().
 */
class Result {

public:

	/*!
	 * @brief Default constructor.
	 */
	Result();

	/*!
	 * @brief Sets the lower bound for GED.
	 * @param[in] lower_bound The lower bound for GED.
	 */
	void set_lower_bound(double lower_bound);

	/*!
	 * @brief Returns the lower bound for GED.
	 * @return The lower bound for GED stored in the result. Equals 0.0 by default.
	 */
	double lower_bound() const;

	/*!
	 * @brief Returns the upper bound for GED.
	 * @return The upper bound for GED stored in the result. Equals std::numeric_limits<double>::infinity() by default.
	 * @note Call sort_node_maps_and_set_upper_bound() to ensure that a call to this method returns the best upper bound.
	 */
	double upper_bound() const;

	/*!
	 * @brief Adds an empty node map to the result.
	 * @param[in] num_nodes_g Number of nodes in first input graph.
	 * @param[in] num_nodes_h Number of nodes in second input graph.
	 * @return The index of the newly added node map.
	 */
	std::size_t add_node_map(std::size_t num_nodes_g, std::size_t num_nodes_h);

	/*!
	 * @brief Adds a node map to the result.
	 * @param[in] node_map The node map to be added to the result.
	 * @return The index of the newly added node map.
	 */
	std::size_t add_node_map(const NodeMap & node_map);

	/*!
	 * @brief Provides access to a node map.
	 * @param[in] index_node_map The index of the node map.
	 * @return Reference to the node map with index @p index_node_map.
	 */
	NodeMap & node_map(std::size_t index_node_map);

	/*!
	 * @brief Checks if a node map is already contained in the vector of all node maps and removes it if this is the case.
	 * @param[in] index_node_map The index of the node map whose non-redundancy should be ensured.
	 * @return Boolean @p true if the node map at @p index_node_map is non-redundant and @p false otherwise.
	 */
	bool is_non_redundant_node_map(std::size_t index_node_map);

	/*!
	 * @brief Provides access to all node maps.
	 * @return Reference to vector of node maps. If sort_node_maps_and_set_upper_bound() has been called, the vector is sorted w.r.t. non-decreasing induced costs.
	 */
	std::vector<NodeMap> & node_maps();

	/*!
	 * @brief Returns the number of node maps.
	 * @return The number of node maps stored in the result.
	 */
	std::size_t num_node_maps() const;

	/*!
	 * @brief Sorts the vector of node maps w.r.t non-decreasing induced cost and possibly discards expensive node maps.
	 * @param[in] num_node_maps The number of node maps to be retained. If smaller than the number of node maps contained in the result, the most expensive node maps are deleted.
	 */
	void sort_node_maps_and_set_upper_bound(std::size_t num_node_maps = std::numeric_limits<std::size_t>::max());

private:

	std::vector<NodeMap> node_maps_;

	double lower_bound_;
};

}

#include "result.ipp"

#endif /* SRC_ENV_RESULT_HPP_ */
