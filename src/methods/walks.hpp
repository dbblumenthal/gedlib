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
 * @file  walks.hpp
 * @brief ged::Walks class declaration.
 */

#ifndef SRC_METHODS_WALKS_HPP_
#define SRC_METHODS_WALKS_HPP_

namespace ged {

/*!
 * @brief Computes an upper bound for general edit costs.
 * @details Implements the method Walks suggested in:
 * - B. Ga&uuml;z&egrave;re, S. Bougleux, K. Riesen, and L. Brun:
 *   &ldquo;Approximate graph edit distance guided by bipartite matching of bags of walks&rdquo;,
 *   https://doi.org/10.1007/978-3-662-44415-3_8
 *
 * Supports the following option in addition to the ones supported by ged::LSAPEBasedMethod:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--load @<filename@></tt> | path to existing configuration file | not specified | n.a. |
 * | <tt>\--save @<filename@></tt> | path where to save configuration file | not specified | n.a. |
 * | <tt>\--depth-range @<smaller convertible to int greater 0@>,@<larger convertible to int greater 0@></tt> | range that specifies possible depths of the walks | <tt>1,5</tt> | if the range is larger than one, the best choice is determined via cross-validation during initialization and the mean is used if the method is run without prior initialization |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Walks : public LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Walks();

	Walks(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	typedef std::map<LabelID, double> Histogram_;

	typedef std::map<LabelID, std::vector<std::size_t>> InverseLabelIndex_;

	typedef std::pair<GEDGraph::NodeID, GEDGraph::NodeID> ProductNode_;

	typedef std::map<std::pair<GEDGraph::NodeID, GEDGraph::NodeID>, std::size_t> ProductNodeSizeTMap_;

	static std::size_t undefined_() {return std::numeric_limits<std::size_t>::max();}

	class AdjGraph_ {
	public:
		AdjGraph_(const GEDGraph & graph);

		AdjGraph_(const AdjGraph_ & adj_graph);

		AdjGraph_();

		void operator=(const AdjGraph_ & adj_graph);

		GEDGraph::NodeID node(std::size_t node_id) const;

		std::pair<std::vector<std::size_t>::const_iterator, std::vector<std::size_t>::const_iterator> nodes_with_label(LabelID label) const;

		std::size_t size() const;

		const Histogram_ & num_walks_from_node_to_labels(std::size_t node_id) const;

		double operator() (std::size_t row, std::size_t col) const;

		void compute_num_walks_(const std::set<LabelID> & node_labels, std::size_t depth);

	private:
		IMatrix adj_matrix_;

		std::vector<Histogram_> num_walks_from_nodes_to_labels_;

		GEDGraph::SizeTNodeMap nodes_;

		InverseLabelIndex_ inverse_label_index_;
	};


	class ProductGraph_ {
	public:
		ProductGraph_(const GEDGraph & g, const GEDGraph & h);

		std::pair<std::vector<std::size_t>::const_iterator, std::vector<std::size_t>::const_iterator> nodes_with_label(LabelID label) const;

		std::size_t node_id(GEDGraph::NodeID node_g, GEDGraph::NodeID node_h) const;

		void compute_num_unmatched_walks_(const std::set<LabelID> & node_labels, const AdjGraph_ & adj_graph_g, const AdjGraph_ & adj_graph_h, std::size_t depth);

		double num_substituted_walks_ending_at_same_label(std::size_t row, std::size_t col) const;

		double num_substituted_walks_ending_at_different_labels(std::size_t row, std::size_t col) const;

		double num_inserted_or_deleted_walks(std::size_t row, std::size_t col) const;

	private:
		IMatrix adj_matrix_;

		DMatrix num_substituted_walks_ending_at_same_label_;

		DMatrix num_substituted_walks_ending_at_different_labels_;

		DMatrix num_inserted_or_deleted_walks_;

		std::vector<ProductNode_> nodes_;

		InverseLabelIndex_ inverse_label_index_;

		ProductNodeSizeTMap_ product_node_to_id_;
	};

	std::size_t depth_;

	std::size_t min_depth_;

	std::size_t max_depth_;

	std::string infile_;

	std::string outfile_;

	std::map<GEDGraph::GraphID, AdjGraph_> adj_graphs_;

	// Member functions inherited from LSAPEBasedMethod.

	virtual void lsape_set_default_options_() final;

	virtual std::string lsape_valid_options_string_() const final;

	virtual bool lsape_parse_option_(const std::string & option, const std::string & arg) final;

	virtual void lsape_init_graph_(const GEDGraph & graph) final;

	virtual void lsape_init_() final;

	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) final;

	virtual void lsape_pre_graph_init_(bool called_at_runtime) final;

	// Helper member functions.

	void init_node_labels_(const GEDGraph & g, const GEDGraph & h, std::set<LabelID> & node_labels) const;

	bool load_config_file_() const;

};

}

#endif /* SRC_METHODS_WALKS_HPP_ */
