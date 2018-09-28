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
 * @file result.ipp
 * @brief ged::Result struct definition.
 */

#ifndef SRC_ENV_RESULT_IPP_
#define SRC_ENV_RESULT_IPP_

namespace ged {

Result ::
Result() :
node_maps_(),
lower_bound_{0.0} {}

void
Result ::
set_lower_bound(double lower_bound) {
	lower_bound_ = lower_bound;
}

double
Result ::
lower_bound() const {
	return lower_bound_;
}

double
Result ::
upper_bound() const {
	if (node_maps_.empty()) {
		return std::numeric_limits<double>::infinity();
	}
	return node_maps_.at(0).induced_cost();
}

std::size_t
Result ::
add_node_map(std::size_t num_nodes_g, std::size_t num_nodes_h) {
	node_maps_.emplace_back(num_nodes_g, num_nodes_h);
	return (node_maps_.size() - 1);
}

std::size_t
Result ::
add_node_map(const NodeMap & node_map) {
	node_maps_.emplace_back(node_map);
	return (node_maps_.size() - 1);
}

NodeMap &
Result ::
node_map(std::size_t index_node_map) {
	return node_maps_.at(index_node_map);
}

bool
Result ::
is_non_redundant_node_map(std::size_t index_node_map) {
	for (std::size_t pos{0}; pos < index_node_map; pos++) {
		if (node_maps_.at(pos) == node_maps_.at(index_node_map)) {
			node_maps_.erase(node_maps_.begin() + index_node_map);
			return false;
		}
	}
	for (std::size_t pos{index_node_map + 1}; pos < node_maps_.size(); pos++) {
		if (node_maps_.at(pos) == node_maps_.at(index_node_map)) {
			node_maps_.erase(node_maps_.begin() + index_node_map);
			return false;
		}
	}
	return true;
}

std::vector<NodeMap> &
Result ::
node_maps() {
	return node_maps_;
}

std::size_t
Result ::
num_node_maps() const {
	return node_maps_.size();
}

void
Result ::
sort_node_maps_and_set_upper_bound(std::size_t num_node_maps) {
	if (node_maps_.empty()) {
		return;
	}
	std::sort(node_maps_.begin(), node_maps_.end());
	if (node_maps_.size() > num_node_maps) {
		node_maps_.erase(node_maps_.begin() + num_node_maps + 1, node_maps_.end());
	}
}

}

#endif /* SRC_ENV_RESULT_IPP_ */
