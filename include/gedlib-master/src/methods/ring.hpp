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
 * @file  ring.hpp
 * @brief ged::Ring class declaration
 */

#ifndef SRC_METHODS_RING_HPP_
#define SRC_METHODS_RING_HPP_

namespace ged {

/*!
 * @brief Computes an upper bound for general edit costs.
 * @details Implements the method %Ring suggested in:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun.
 *   &ldquo;Ring-based approximation of graph edit distance&rdquo;,
 *   https://doi.org/10.1007/978-3-319-97785-0_28
 *
 * Supports the following option in addition to the ones supported by ged::LSAPEBasedMethod:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--load @<filename@></tt> | path to existing configuration file | not specified | n.a. |
 * | <tt>\--save @<filename@></tt> | path where to save configuration file | not specified | n.a. |
 * | <tt>\--led-method LSAPE_OPTIMAL\|LSAPE_GREEDY\|GAMMA</tt> | method for computing the layer distances | @p LSAPE_OPTIMAL | see S+SSPR paper <br> if not @p GAMMA, the method @p \--sort-method has no effect |
 * | <tt>\--sort-method STD\|COUNTING</tt> | the employed sorting algorithm | @p COUNTING | @ref ged::util::counting_sort() <br> use counting sort if the number of different edge labels is constant |
 * | <tt>\--init-evaluations @<convertible to int greater 0@></tt> | number of blackbox evaluations during initialization | @p 50 | see S+SSPR paper |
 * | <tt>\--init-num-initial-solutions @<convertible to int greater 0@></tt> | number of initial solutions used for initialization | @p 50 | see S+SSPR paper |
 * | <tt>\--init-mu @<convertible to double between 0 and 1@></tt> | parameter @f$\mu@f$ used for initialization | @p 1 | see S+SSPR paper |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Ring : public LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Ring();

	Ring(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	enum LEDMethod_ {LSAPE_OPTIMAL, LSAPE_GREEDY, GAMMA};

	enum SortMethod_ {STD, COUNTING};

	struct Layer_ {
		Layer_(std::size_t level);

		std::size_t level;

		std::vector<LabelID> node_labels;

		std::vector<LabelID> inner_edge_labels;

		std::vector<LabelID> outer_edge_labels;
	};

	struct Ring_ {
		Ring_();

		std::vector<Layer_> layers;
	};

	class Evaluator_ : public NOMAD::Evaluator {
	public:
		Evaluator_(const NOMAD::Parameters & param, Ring<UserNodeLabel, UserEdgeLabel> * ring);

		~Evaluator_();

		bool eval_x(NOMAD::Eval_Point & x, const NOMAD::Double & h_max, bool & count_eval) const;

	private:
		Ring<UserNodeLabel, UserEdgeLabel> * ring_;
	};

	typedef std::map<GEDGraph::NodeID, Ring_> NodeRingMap_;

	std::map<GEDGraph::GraphID, NodeRingMap_> rings_;

	LEDMethod_ led_method_;

	SortMethod_ sort_method_;

	std::size_t num_layers_;

	std::vector<double> alpha_;

	std::vector<double> lambda_;

	double mu_;

	std::size_t num_evals_;

	std::size_t num_x0s_;

	std::string infile_;

	std::string outfile_;

	// Member functions inherited from LSAPEBasedMethod.

	virtual void lsape_set_default_options_() final;

	virtual std::string lsape_valid_options_string_() const final;

	virtual bool lsape_parse_option_(const std::string & option, const std::string & arg) final;

	virtual void lsape_init_() final;

	virtual void lsape_init_graph_(const GEDGraph & graph) final;

	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) final;

	virtual void lsape_default_post_graph_init_() final;

	virtual void lsape_pre_graph_init_(bool called_at_runtime) final;

	// Private helper member functions.

	void populate_instance_with_params_(const GEDGraph & g, const GEDGraph & h, const vector<double> & alpha, const vector<double> & lambda, DMatrix & lsape_instance) const;

	void set_num_layers_();

	void build_rings_(const GEDGraph & graph);

	void build_ring_(const GEDGraph & graph, GEDGraph::NodeID root, NodeRingMap_ & rings);

	void init_x0s_(std::vector<NOMAD::Point> & x0s) const;

	void nomad_point_to_params_(const NOMAD::Point & x, std::vector<double> & alpha, std::vector<double> & lambda) const;

	void normalize_params_();

	void eval_x_(NOMAD::Eval_Point & x) const;

	bool load_config_file_() const;

	double compute_ring_distance_(const GEDGraph & g, const GEDGraph & h, const NodeRingMap_ & rings_g, const NodeRingMap_ & rings_h,
				const std::vector<double> & alpha, const std::vector<double> & lambda, std::size_t row_in_master, std::size_t col_in_master) const;

	double compute_substitution_cost_(const Ring_ & ring_i, const Ring_ & ring_k, const std::vector<double> & alpha, std::size_t level) const;

	double compute_deletion_cost_(const Ring_ & ring, const std::vector<double> & alpha, std::size_t level) const;

	double compute_insertion_cost_(const Ring_ & ring, const std::vector<double> & alpha, std::size_t level) const;

	double compute_layer_distance_(const Layer_ & lhs, const Layer_ & rhs, const std::vector<double> & alpha) const;

	double lsape_multiset_cost_(const std::vector<LabelID> & lhs, const std::vector<LabelID> & rhs, bool node_labels) const;

	double gamma_multiset_cost_(const std::vector<LabelID> & lhs, const std::vector<LabelID> & rhs, bool node_labels) const;

	void write_params_to_file_() const;

	void read_params_from_file_();
};

}

#endif /* SRC_METHODS_RING_HPP_ */
