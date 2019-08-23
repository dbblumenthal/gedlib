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
 * @file  ring_ml.hpp
 * @brief ged::RingML class declaration.
 */

#ifndef SRC_METHPODS_RING_ML_HPP_
#define SRC_METHPODS_RING_ML_HPP_

namespace ged {

/*!
 * @brief Uses ring structures for defining feature vectors for node edit operations.
 * @details Implements the ring based feature vectors suggested in:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun.
 *   &ldquo;Upper bounding GED via transformations to LSAPE based on rings and machine learning.&rdquo;,
 *   to be submitted to TKDE
 *
 * Supports the following option in addition to the ones supported by ged::MLBasedMethod:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--led-method LSAPE_OPTIMAL\|LSAPE_GREEDY\|GAMMA</tt> | method for computing the layer distances | @p LSAPE_OPTIMAL | see PR paper <br> if not @p GAMMA, the method @p \--sort-method has no effect |
 * | <tt>\--sort-method STD\|COUNTING</tt> | the employed sorting algorithm | @p COUNTING | @ref ged::util::counting_sort() <br> use counting sort if the number of different edge labels is constant |
 * | <tt>\--topological-features TRUE\|FALSE</tt> | decides whether to include features for representing the topologies of the graphs | @p TRUE | see PR paper |
 * | <tt>\--global-features TRUE\|FALSE</tt> | decides whether to include features for representing global properties of the graphs | @p TRUE | see PR paper |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class RingML : public MLBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~RingML();

	RingML(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

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

	typedef std::map<GEDGraph::NodeID, Ring_> NodeRingMap_;

	std::map<GEDGraph::GraphID, NodeRingMap_> rings_;

	LEDMethod_ led_method_;

	SortMethod_ sort_method_;

	bool use_topological_features_;

	bool use_global_features_;

	std::size_t num_layers_;

	std::vector<double> global_features_;

	// Member functions inherited from MLBasedMethod.

	virtual void ml_init_graph_(const GEDGraph & graph) final;

	virtual void ml_set_default_options_() final;

	virtual std::string ml_valid_options_string_() const final;

	virtual bool ml_parse_option_(const std::string & option, const std::string & arg) final;

	virtual void ml_init_feature_variables_(const GEDGraph & g, const GEDGraph & h, std::size_t num_threads) final;

	virtual void ml_populate_substitution_feature_vector_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k, std::vector<double> & feature_vector) final;

	virtual void ml_populate_deletion_feature_vector_(const GEDGraph & g, GEDGraph::NodeID i, std::vector<double> & feature_vector) final;

	virtual void ml_populate_insertion_feature_vector_(const GEDGraph & h, GEDGraph::NodeID k, std::vector<double> & feature_vector) final;

	virtual std::size_t ml_get_num_features_() final;

	virtual void ml_init_for_num_features_() final;

	// Private helper member functions.

	void set_num_layers_();

	void build_rings_(const GEDGraph & graph);

	void build_ring_(const GEDGraph & graph, GEDGraph::NodeID root, NodeRingMap_ & rings);

	void add_global_features_(std::vector<double> & feature_vector) const;

	void add_layer_substitution_features_(const Ring_ & ring_i, const Ring_ & ring_k, std::size_t level, std::vector<double> & feature_vector) const;

	void add_layer_deletion_features_(const Ring_ & ring, std::size_t level, std::vector<double> & feature_vector) const;

	void add_layer_insertion_features_(const Ring_ & ring, std::size_t level, std::vector<double> & feature_vector) const;

	void add_layer_features_(const Layer_ & lhs, const Layer_ & rhs, std::vector<double> & feature_vector) const;

	double lsape_multiset_cost_(const std::vector<LabelID> & lhs, const std::vector<LabelID> & rhs, bool node_labels) const;

	double gamma_multiset_cost_(const std::vector<LabelID> & lhs, const std::vector<LabelID> & rhs, bool node_labels) const;
};

}

#endif /* SRC_METHPODS_RING_ML_HPP_ */
