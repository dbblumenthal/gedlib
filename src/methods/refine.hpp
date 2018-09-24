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
 * Does not support any options except for the ones supported by ged::LSBasedMethod.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Refine : public LSBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Refine();

	Refine(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	virtual void ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) final;

};

}

#endif /* SRC_METHODS_REFINE_HPP_ */
