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
 * @file bipartite_ml.hpp
 * @brief ged::BipartiteML class declaration.
 */

#ifndef SRC_METHODS_BIPARTITE_ML_HPP_
#define SRC_METHODS_BIPARTITE_ML_HPP_

namespace ged {

/*!
 * @brief Uses characteristics of an LSAPE instance for defining feature vectors for node edit operations.
 * @details Implements the feature vectors suggested in:
 * - K. Riesen and M. Ferrer:
 *   &ldquo;Predicting the correctness of node assignments in bipartite graph matching&rdquo;,
 *   https://doi.org/10.1016/j.patrec.2015.10.007,
 *
 * and their extension to general instantiations of LSAPE-GED suggested in:
 * - D. B. Blumenthal, S. Bougleux, J. Gamper, L. Brun:
 *   &ldquo;Upper bounding GED via transformations to LSAPE based on rings and machine learning.&rdquo;,
 *   To be submitted to TKDE.
 *
 * Supports the following options in addition to the ones supported by ged::MLBasedMethod.
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--lsape-method BIPARTITE\|BRANCH_FAST\|BRANCH_UNIFORM\|BRANCH\|NODE\|RING\|SUBGRAPH\|WALKS</tt> | method for populating the LSAPE instance | @p BIPARTITE | if @p BIPARTITE, the feature vectors are identical to the ones suggested in https://doi.org/10.1016/j.patrec.2015.10.007 |
 * | <tt>\--lsape-options '[--@<option@> @<arg@>] [...]'</tt> | options string passed to the method used for populating the LSAPE instance | @p '' | ged::Bipartite, ged::BranchFast, ged::BranchUniform, ged::Branch, ged::Node, ged::Ring, ged::Subgraph, ged::Walks |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class BipartiteML : public MLBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~BipartiteML();

	BipartiteML(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	typedef Eigen::Array<double, 1, Eigen::Dynamic> RowVector_;

	typedef Eigen::Array<double, Eigen::Dynamic, 1> ColumnVector_;

	class RowFeatures_ {

	public:

		RowFeatures_();

		void init(const Eigen::ArrayXXd & substitution_matrix);

		void add_features_(const Eigen::ArrayXXd & matrix, std::size_t row, std::size_t col, std::vector<double> & feature_vector) const;

	private:

		RowVector_ maxima_;

		RowVector_ minima_;

		RowVector_ means_;

		RowVector_ deviations_;

		RowVector_ leaders_;

		RowVector_ intervals_;
	};

	class ColFeatures_ {

	public:

		ColFeatures_();

		void init(const Eigen::ArrayXXd & substitution_matrix);

		void add_features_(const Eigen::ArrayXXd & matrix, std::size_t row, std::size_t col, std::vector<double> & feature_vector) const;

	private:

		ColumnVector_ maxima_;

		ColumnVector_ minima_;

		ColumnVector_ means_;

		ColumnVector_ deviations_;

		ColumnVector_ leaders_;

		ColumnVector_ intervals_;
	};

	LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> * lsape_method_;

	std::string lsape_method_options_;

	Eigen::ArrayXXd lsape_instance_;

	std::vector<double> global_features_;

	RowFeatures_ row_features_;

	ColFeatures_ col_features_;

	// Member functions inherited from MLBasedMethod.

	virtual void ml_init_() final;

	virtual void ml_set_default_options_() final;

	virtual bool ml_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ml_valid_options_string_() const final;

	virtual void ml_init_feature_variables_(const GEDGraph & g, const GEDGraph & h, std::size_t num_threads) final;

	virtual void ml_populate_substitution_feature_vector_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k, std::vector<double> & feature_vector) final;

	virtual void ml_populate_deletion_feature_vector_(const GEDGraph & g, GEDGraph::NodeID i, std::vector<double> & feature_vector) final;

	virtual void ml_populate_insertion_feature_vector_(const GEDGraph & h, GEDGraph::NodeID k, std::vector<double> & feature_vector) final;

	virtual std::size_t ml_get_num_features_() final;

	// Private helper member functions.

	void populate_lsape_instance_(const GEDGraph & g, const GEDGraph & h, std::size_t num_threads);

	void compute_global_features_(const Eigen::ArrayXXd & substitution_matrix);

	void add_global_features_(std::vector<double> & feature_vector) const;

	void add_cell_features_(std::size_t row, std::size_t col, double node_cost, std::vector<double> & feature_vector) const;

};

}

#endif /* SRC_METHODS_BIPARTITE_ML_HPP_ */
