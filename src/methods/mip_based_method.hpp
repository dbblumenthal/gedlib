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
 * @file mip_based_method.hpp
 * @brief ged::MIPBasedMethod class declaration.
 */

#ifndef SRC_METHODS_MIP_BASED_METHOD_HPP_
#define SRC_METHODS_MIP_BASED_METHOD_HPP_

namespace ged {

/*!
 * @brief Abstract class for methods that use mixed integer linear programming for exactly or approximatively computing the graph edit distance.
 * @details All derived methods support the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--threads @<convertible to int greater 0@></tt> | number of threads | @p 1 | number of threads to be used by MIP/LP solver |
 * | <tt>\--time-limit @<convertible to double@></tt> | time limit in seconds | @p 0 | if less or equal @p 0, no time limit is enforced |
 * | <tt>\--relax TRUE\|FALSE </tt> | if @p TRUE, all integrality constraints are relaxed | @p FALSE | if @p TRUE, the model populated by mip_populate_model_() must be continuous |
 * | <tt>\--project-to-node-map TRUE\|FALSE </tt> | if @p TRUE, continuous solutions of relaxed models are projected to node maps | @p TRUE | mip_model_to_lsape_projection_problem_() |
 * | <tt>\--tune TRUE\|FALSE </tt> | if @p TRUE, the parameters of the model are tuned before optimization | @p FALSE | <a href="http://www.gurobi.com/documentation/8.0/refman/cpp_grbmodel_tune.html#cppmethod:GRBModel::tune">Gurobi Manual</a> |
 * | <tt>\--tune-time-limit @<convertible to double@></tt> | time limit for parameter tuning in seconds | @p 0 | if less or equal @p 0, the time limit is set automatically |
 * | <tt>\--map-root-to-root TRUE\|FALSE</tt> | decide if the roots of the input graphs are mapped to each other | @p FALSE | if @p TRUE, the nodes with ID 0 of the input graphs are mapped to each other |
 * | <tt>\--lsape-model ECBP\|EBP\|FLWC\|FLCC\|FBP\|SFBP\|FBP0</tt> | model for optimally solving LSAPE | @p ECBP | ged::LSAPESolver::Model |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class MIPBasedMethod : public GEDMethod<UserNodeLabel, UserEdgeLabel> {

public:

	/*!
	 * @brief Pure virtual destructor.
	 * @note Must be implemented by derived classes.
	 */
	virtual ~MIPBasedMethod() = 0;

	/*!
	 * @brief Constructor.
	 * @param[in] ged_data The instance on which the method should be run.
	 */
	MIPBasedMethod(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

protected:

	/*!
	 * @brief If @p true, the model populated by mip_populate_model_() must be continuous.
	 */
	bool relax_;

	/*!
	 * @brief If @p true, the model populated by mip_populate_model_() must enforce that the nodes with ID 0 are mapped to each other.
	 */
	bool map_root_to_root_;

private:

	std::size_t num_threads_;

	double time_limit_in_sec_;

	bool tune_;

	double tune_time_limit_in_sec_;

	LSAPESolver::Model lsape_model_;

	bool project_to_node_map_;

	// Member functions inherited from GEDMethod.

	virtual void ged_init_() final;

	virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;

	virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string ged_valid_options_string_() const final;

	virtual void ged_set_default_options_() final;

	// Virtual member functions to be overridden by derived classes.

	/*!
	 * @brief Runs the local search from an initial node map.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[out] model The mixed integer linear programming model. Must be continuous if @p relax_ equals @p true.
	 * @note Must be overridden by derived classes of ged::MIPBasedMethod.
	 */
	virtual void mip_populate_model_(const GEDGraph & g, const GEDGraph & h, GRBModel & model);

	/*!
	 * @brief Given a, possibly sub-optimally, solved unrelaxed model, this method constructs a node map and sets its induced cost.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[in] model Possibly sub-optimally solved unrelaxed model.
	 * @param[out] node_map Node map which has to be constructed and whose induced cost has to be set.
	 * @note Must be overridden by derived classes of ged::MIPBasedMethod.
	 */
	virtual void mip_model_to_node_map_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, NodeMap & node_map);

	/*!
	 * @brief Given a, possibly sub-optimally, solved model, this method constructs an LSAPE instance for projecting a continuous solution to a node map.
	 * @param[in] g Input graph.
	 * @param[in] h Input graph.
	 * @param[in] model Possibly sub-optimally solved model.
	 * @param[in,out] lsape_instance The LSAPE instance to be constructed. Initialized as matrix of ones.
	 * @return Boolean @p false by default. Derived classes that override this method must return @p true.
	 * @note Must be overridden by derived classes of ged::MIPBasedMethod that want to provide upper bounds when called with the option <tt>--relax TRUE</tt>.
	 * The LSAPE instance @p lsape_instance should be constructed such that, when constructed for an unrelaxed model, the unique optimal solution for
	 * @p lsape_instance corresponds to the node map constructed by mip_model_to_node_map_().
	 */
	virtual bool mip_model_to_lsape_projection_problem_(const GEDGraph & g, const GEDGraph & h, GRBModel & model, DMatrix & lsape_instance);

	/*!
	 * @brief Initializes the method.
	 * @note Must be overridden by derived classes of ged::MIPBasedMethod that require initialization.
	 */
	virtual void mip_init_();

	/*!
	 * @brief Parses one option that is not among the ones shared by all derived classes of ged::MIPBasedMethod.
	 * @param[in] option The name of the option.
	 * @param[in] arg The argument of the option.
	 * @return Returns true if @p option is a valid option name for the method and false otherwise.
	 * @note Must be overridden by derived classes of ged::MIPBasedMethod that have options that are not among the ones shared by all derived classes of ged::MIPBasedMethod.
	 */
	virtual bool mip_parse_option_(const std::string & option, const std::string & arg);

	/*!
	 * @brief Returns string of all valid options that are not among the ones shared by all derived classes of ged::MIPBasedMethod.
	 * @return String of the form @"[--@<option@> @<arg@>] [...]@".
	 * @note Must be overridden by derived classes of ged::MIPBasedMethod that have options that are not among the ones shared by all derived classes of ged::MIPBasedMethod.
	 */
	virtual std::string mip_valid_options_string_() const;

	/*!
	 * @brief Sets all options that are not among the ones shared by all derived classes of ged::MIPBasedMethod to default values.
	 * @note Must be overridden by derived classes of ged::MIPBasedMethod that have options that are not among the ones shared by all derived classes of ged::MIPBasedMethod.
	 */
	virtual void mip_set_default_options_();

};

}

#endif /* SRC_METHODS_MIP_BASED_METHOD_HPP_ */
