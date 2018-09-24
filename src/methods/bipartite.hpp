/*!
 * @file  bipartite.hpp
 * @brief ged::Bipartite class declaration.
 */

#ifndef SRC_METHODS_BIPARTITE_HPP_
#define SRC_METHODS_BIPARTITE_HPP_

namespace ged {

/*!
 * @brief Computes lower and upper bounds for general edit costs.
 * @details Implements the method Bipartite suggested in:
 * - K. Riesen and H. Bunke:
 *   &ldquo;Approximate graph edit distance computation by means of bipartite graph matching&rdquo;,
 *   https://doi.org/10.1016/j.imavis.2008.04.004
 *
 * Does not support any options except for the ones supported by ged::LSAPEBasedMethod.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Bipartite : public LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Bipartite();

	Bipartite(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	// Inherited member functions from LSAPEBasedMethod.

	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance) final;

	// Helper member functions.

	double compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k) const;

	double compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const;

	double compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const;
};

}

#endif /* SRC_METHODS_BIPARTITE_HPP_ */
