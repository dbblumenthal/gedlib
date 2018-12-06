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
 * @file refine.hpp
 * @brief ged::Refine class declaration.
 */

#ifndef SRC_METHODS_REFINE_HPP_
#define SRC_METHODS_REFINE_HPP_

namespace ged {

/*!
 * @brief Computes an upper bound for general edit costs.
 * @details Implements the method Refine suggested in:
 * - Z. Zeng, A. K. H. Tung, J. Wang, J. Feng, and L. Zhou:
 *   &ldquo;Comparing stars: On approximating graph edit distance&rdquo;,
 *   http://dx.doi.org/10.14778/1687627.1687631
 *
 * Supports the following options in addition tothe ones supported by ged::LSBasedMethod.
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--max-swap-size @<convertible to int greater equal 2@> | maximum size of swap | @p 2 | If set to a value greater @p 2, K-Refine is used. |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Refine : public LSBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Refine();

	Refine(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	std::size_t max_swap_size_;

	struct Swap_{
		std::vector<NodeMap::Assignment> original_assignments;
		std::vector<NodeMap::Assignment> new_assignments;

		double cost(const GEDGraph & g, const GEDGraph & h, const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data, NodeMap & node_map) const;

		void do_swap(NodeMap & node_map, double delta_cost=0) const;

		void undo_swap(NodeMap & node_map) const;

		std::string print() const;

	};

	virtual void ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map,NodeMap & output_node_map) final;

	virtual bool ls_parse_option_(const std::string & option, const std::string & arg);

	virtual std::string ls_valid_options_string_() const;

	bool next_subset_(const std::size_t set_size, std::vector<std::size_t> & subset);

	virtual void ls_set_default_options_();

};

}

#endif /* SRC_METHODS_REFINE_HPP_ */
