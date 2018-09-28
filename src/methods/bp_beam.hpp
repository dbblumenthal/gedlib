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
 * @file bp_beam.hpp
 * @brief ged::BPBeam class declaration.
 */

#ifndef SRC_METHODS_BP_BEAM_HPP_
#define SRC_METHODS_BP_BEAM_HPP_

namespace ged {

/*!
 * @brief Computes an upper bound for general edit costs.
 * @details Implements the methods BP-Beam suggested in:
 * - K. Riesen, A. Fischer, and H. Bunke:
 *   &ldquo;Combining bipartite graph matching and beam search for graph edit distance approximation&rdquo;,
 *   https://doi.org/10.1007/978-3-319-11656-3_11,
 *
 * the extension IBP-Beam suggested in:
 * - M. Ferrer, F. Serratosa, and K. Riesen:
 *   &ldquo;A first step towards exact graph edit distance using bipartite graph matching&rdquo;,
 *   https://doi.org/10.1007/978-3-319-18224-7_8,
 *
 * and the extension MLBP-Beam suggested in:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun:
 *   &ldquo;Designing heuristics for graph edit distance.&rdquo;,
 *   To be submitted to PR.
 *
 * Supports the following options in addition to the ones supported by ged::LSBasedMethod:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * |<tt>\--beam-size @<convertible to int greater 0@></tt> | the size of the priority queue used by the beam search | @p 5 | https://doi.org/10.1007/978-3-319-11656-3_11 |
 * |<tt>\--ordering-method BIPARTITE_ML\|RING_ML\|RANDOM</tt> | the machine learning method used for ordering the assignments | @p RANDOM | if @p BIPARTITE_ML or @p RING_ML, MLBP-BEAM is used and the option @p \--num-random-orderings has no effect <br> if @p RANDOM, the option @p \--ordering-options has no effect |
 * |<tt>\--ordering-options '[--@<option@> @<arg@>] [...]'</tt> | the options string passed to the machine learning method | @p '' | ged::BipartiteML and ged::RingML |
 * |<tt>\--num-random-orderings @<convertible to int greater equal 0@></tt> | the number of random orderings | @p 0 | https://doi.org/10.1007/978-3-319-18224-7_8 <br> if greater 0 and not overruled by options @p \--ordering-method, IBP-Beam is used |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class BPBeam : public LSBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~BPBeam();

	BPBeam(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	class TreeNode_ {

	public:

		TreeNode_(const TreeNode_ & tree_node);

		TreeNode_(const NodeMap & node_map);

		TreeNode_(const TreeNode_ & tree_node, const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data, const GEDGraph & g, const GEDGraph & h, const NodeMap::Assignment & assignment_t, const NodeMap::Assignment & assignment_2);

		const NodeMap & node_map() const;

		std::size_t depth() const;

		bool operator<(const TreeNode_ & tree_node) const;

	private:

		NodeMap node_map_;

		std::size_t depth_;
	};

	std::size_t beam_size_;

	std::size_t num_random_orderings_;

	MLBasedMethod<UserNodeLabel, UserEdgeLabel> * ordering_method_;

	// Member functions inherited from LSBasedMethod.

	virtual void ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) final;

	virtual void ls_set_default_options_() final;

	virtual bool ls_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ls_valid_options_string_() const final;

	virtual void ls_init_() final;

	// Private helper member functions.

	void sort_and_shrink_to_beam_size_(std::vector<TreeNode_> & open) const;

	void order_assignments_(const GEDGraph & g, const GEDGraph & h, std::vector<NodeMap::Assignment> & assignments);

};

}




#endif /* SRC_METHODS_BP_BEAM_HPP_ */
