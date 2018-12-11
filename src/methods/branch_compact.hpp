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
 *  @file  branch_compact.hpp
 *  @brief ged::BranchCompact class declaration.
 */

#ifndef SRC_METHODS_BRANCH_COMPACT_HPP_
#define SRC_METHODS_BRANCH_COMPACT_HPP_

namespace ged {

/*!
 * @brief Computes a lower bound for uniform edit costs.
 * @details Implements the method %BranchCompact suggested in:
 * - W. Zheng, L. Zou, X. Lian, D. Wang, and D. Zhao:
 *   &ldquo;Efficient graph similarity search over large graph databases&rdquo;,
 *   https:://doi.org/10.1109/TKDE.2014.2349924
 *
 * Supports the following option:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--sort-method STD\|COUNTING</tt> | the employed sorting algorithm | @p COUNTING | @ref ged::util::counting_sort() <br> use counting sort if the number of different edge labels is constant |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class BranchCompact : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~BranchCompact();

	BranchCompact(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	enum SortMethod_ {STD, COUNTING};

	class SortedUserEdgeLabels_ {
	public:
		SortedUserEdgeLabels_(const GEDGraph & g, SortMethod_ sort_method);

		const std::vector<LabelID> & get_incident_labels(GEDGraph::NodeID) const;

	private:
		std::map<GEDGraph::NodeID, std::vector<LabelID>> sorted_edge_labels_;
	};

	struct Branch_ {

		Branch_(LabelID node_label, const std::vector<LabelID> & sorted_edge_labels);

		Branch_(const Branch_ & branch);

		LabelID node_label;

		std::vector<LabelID> sorted_edge_labels;

		int compare(const Branch_ & rhs) const;

		bool operator<(const Branch_ & rhs) const;

		bool operator>(const Branch_ & rhs) const;

		bool operator==(const Branch_ & rhs) const;
	};

	SortMethod_ sort_method_;

	std::map<GEDGraph::GraphID, std::list<Branch_>> branches_;

	// Member functions inherited from GEDMethod.

	virtual void ged_init_() final;

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ged_valid_options_string_() const final;

	virtual void ged_set_default_options_() final;

	// Private helper member functions.

	void init_graph_(const GEDGraph & graph);
};

}

#endif /* SRC_METHODS_BRANCH_COMPACT_HPP_ */
