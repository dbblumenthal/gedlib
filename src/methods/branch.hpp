/*!
 * @file  branch.hpp
 * @brief ged::Branch class declaration.
 */

#ifndef SRC_METHODS_BRANCH_HPP_
#define SRC_METHODS_BRANCH_HPP_

namespace ged {

/*!
 * @brief Computes lower and upper bounds for general edit costs.
 * @details Implements the method Branch suggested in:
 * - D. B. Blumenthal and J. Gamper:
 *   &ldquo;Improved lower bounds for graph edit distance&rdquo;,
 *   https:://doi.org/10.1109/TKDE.2017.2772243
 *
 * Does not support any options except for the ones supported by ged::LSAPEBasedMethod.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Branch : public LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Branch();

	Branch(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	// Inherited member functions from LSAPEBasedMethod.

	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & master_problem) final;

	// Helper member functions.

	double compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const;

	double compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const;

	double compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const;
};

}

#endif /* SRC_METHODS_BRANCH_HPP_ */
