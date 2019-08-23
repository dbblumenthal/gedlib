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
 * @file  subgraph.hpp
 * @brief ged::Subgraph class declaration.
 */

#ifndef SRC_METHODS_SUBGRAPH_HPP_
#define SRC_METHODS_SUBGRAPH_HPP_

namespace ged {

/*!
 * @brief Computes upper bounds for general edit costs.
 * @details Implements the method %Subgraph suggested in:
 * - V. Carletti, B. Ga&uuml;z&egrave;re, L. Brun, and M. Vento:
 *   &ldquo;Approximate graph edit distance computation combining bipartite matching and exact neighborhood substructure distance&rdquo;,
 *   https://doi.org/10.1007/978-3-319-18224-7_19
 *
 * Supports the following option in addition to the ones supported by ged::LSAPEBasedMethod:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--load @<filename@></tt> | path to existing configuration file | not specified | n.a. |
 * | <tt>\--save @<filename@></tt> | path where to save configuration file | not specified | n.a. |
 * | <tt>\--subproblem-solver ANCHOR_AWARE_GED\|%F1\|%F2\|COMPACT_MIP</tt> | method for exactly solving the subproblems | @p ANCHOR_AWARE_GED | the methods %F1, %F2, and COMPACT_MIP are available only if GEDLIB is installed with Gurobi |
 * | <tt>\--subproblem-solver-options '[--@<option@> @<arg@>] [...]'</tt> | options string passed to the ground truth method | @p '' | ged::AnchorAwareGED, ged::F1, ged::F2, ged::CompactMIP |
 * | <tt>\--depth-range @<smaller convertible to int greater 0@>,@<larger convertible to int greater 0@></tt> | range that specifies possible depths of the subgraphs | <tt>1,5</tt> | if the range is larger than one, the best choice is determined via cross-validation during initialization and the mean is used if the method is run without prior initialization |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Subgraph : public LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Subgraph();

	Subgraph(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	std::size_t depth_;

	std::size_t min_depth_;

	std::size_t max_depth_;

	std::string infile_;

	std::string outfile_;

	std::string exact_options_;

	std::map<GEDGraph::GraphID, GEDGraph> subgraphs_;

	Options::GEDMethod exact_method_;

	// Inherited member functions from LSAPEBasedMethod.

	virtual void lsape_set_default_options_() final;

	virtual std::string lsape_valid_options_string_() const final;

	virtual bool lsape_parse_option_(const std::string & option, const std::string & arg) final;

	virtual void lsape_init_graph_(const GEDGraph & graph) final;

	virtual void lsape_init_() final;

	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) final;

	virtual void lsape_pre_graph_init_(bool called_at_runtime) final;

	// Helper member functions.

	double compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k, GEDMethod<UserNodeLabel, UserEdgeLabel> * exact_method) const;

	double compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const;

	double compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const;

	void build_subgraphs_(const GEDGraph & graph);

	void build_subgraph_(const GEDGraph & graph, GEDGraph::NodeID node);

	GEDGraph::GraphID subgraph_id_(const GEDGraph & graph, GEDGraph::NodeID node) const;

	void build_subgraph_dfs_(const GEDGraph & graph, GEDGraph::NodeID current_node, std::size_t depth_current_node, GEDGraph::NodeNodeMap & ids_in_subgraph, GEDGraph & subgraph);

	bool load_config_file_() const;

	GEDMethod<UserNodeLabel, UserEdgeLabel> * exact_method_factory_() const;

};

}

#endif /* SRC_METHODS_SUBGRAPH_HPP_ */
