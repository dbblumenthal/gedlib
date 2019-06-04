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
 * @file  ipfp.hpp
 * @brief ged::IPFP class declaration.
 */

#ifndef SRC_METHODS_IPFP_HPP_
#define SRC_METHODS_IPFP_HPP_

namespace ged {

/*!
 * @brief Computes an upper bound for general edit costs.
 * @details Implements the method %IPFP (Frank-Wolfe for QAP and GED) suggested in:
 * - S. Bougleux, L. Brun, V. Carletti, P. Foggia, B. Ga&uuml;z&egrave;re, and M. Vento:
 *   &ldquo;Graph edit distance as a quadratic assignment problem&rdquo;
 *   https://doi.org/10.1016/j.patrec.2016.10.001,
 *
 * the reduction QAPE suggested in:
 * - S. Bougleux, B. Ga&uuml;z&egrave;re, and L. Brun:
 *   &ldquo;Graph edit distance as a quadratic program&rdquo;
 *   https://doi.org/10.1109/ICPR.2016.7899881,
 *
 * the extension m-IPFP suggested in:
 * - &Eacute;. Daller, S. Bougleux, B. Ga&uuml;z&egrave;re, and L. Brun:
 *   &ldquo;Approximate graph edit distance by several local searches in parallel&rdquo;,
 *   https://doi.org/10.5220/0006599901490158,
 *
 * and the extension C-QAP suggested in:
 * - D. B. Blumenthal, &Eacute;. Daller, S. Bougleux, and L. Brun:
 *   &ldquo;Quasimetric graph edit distance as a compact quadratic assignment problem&rdquo;,
 *   accepted for publication at ICPR 2018
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--iterations @<convertible to int greater equal 0@></tt> | maximal number of iterations | @p 100 | if @p 0, no iteration based termination criterion is used |
 * | <tt>\--time-limit @<convertible to double@></tt> | time limit in seconds | @p 0 | if less or equal @p 0, no time limit is enforced |
 * | <tt>\--epsilon @<convertible to double@></tt> | convergence threshold | @p 0.001 | if less or equal @p 0, no convergence based termination criterion is used |
 * | <tt>\--lsape-model ECBP\|EBP\|FLWC\|FLCC\|FBP\|SFBP\|FBP0</tt> | model for optimally solving LSAPE | @p ECBP | ged::LSAPESolver::Model |
 * | <tt>\--quadratic-model QAPE\|B-QAP\|C-QAP</tt> | quadratic assignment model | @p QAPE | ICPR paper |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class IPFP : public LSBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~IPFP();

	IPFP(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

private:

	enum QuadraticModel_ {C_QAP, B_QAP, QAPE, QAP};

	class QAPInstance_ {

	public:

		QAPInstance_();

		void init(const IPFP<UserNodeLabel, UserEdgeLabel> * ipfp, const GEDGraph & g, const GEDGraph & h);

		double operator() (std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const;

		double operator() (std::size_t row, std::size_t col) const;

		std::size_t num_rows() const;

		std::size_t num_cols() const;

		std::size_t num_nodes_g() const;

		std::size_t num_nodes_h() const;

	private:

		double quadratic_cost_qap_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const;

		double quadratic_cost_b_qap_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const;

		double quadratic_cost_qape_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const;

		double quadratic_cost_c_qap_(std::size_t row_1, std::size_t col_1, std::size_t row_2, std::size_t col_2) const;

		const IPFP<UserNodeLabel, UserEdgeLabel> * ipfp_;

		const GEDGraph * g_;

		const GEDGraph * h_;

		std::size_t num_nodes_g_;

		std::size_t num_nodes_h_;

		double translation_factor_;
	};

	QuadraticModel_ quadratic_model_;

	LSAPESolver::Model lsape_model_;

	double epsilon_;

	std::size_t max_itrs_;

	double time_limit_in_sec_;

	double omega_;

	QAPInstance_ qap_instance_;

	// Member functions inherited from LSBasedMethod.

	void ls_run_from_initial_solution_(const GEDGraph & g, const GEDGraph & h, double lower_bound, const NodeMap & initial_node_map, NodeMap & output_node_map) final;

	virtual void ls_runtime_init_(const GEDGraph & g, const GEDGraph & h) final;

	virtual void ls_set_default_options_() final;

	virtual bool ls_parse_option_(const std::string & options, const std::string & arg) final;

	virtual std::string ls_valid_options_string_() const final;

	// Private helper member functions.

	void node_map_to_matrix_(const NodeMap & node_map, DMatrix & matrix) const;

	double compute_induced_linear_cost_(const QAPInstance_ & qap_instance, const DMatrix & x) const;

	double compute_induced_linear_cost_(const QAPInstance_ & qap_instance, const LSAPSolver & solver) const;

	double compute_induced_linear_cost_(const QAPInstance_ & qap_instance, const LSAPESolver & solver) const;

	double compute_induced_quadratic_cost_(const QAPInstance_ & qap_instance, const LSAPSolver & solver) const;

	double compute_induced_quadratic_cost_(const QAPInstance_ & qap_instance, const LSAPESolver & solver) const;

	void solve_linear_problem_(const QAPInstance_ & qap_instance, LSAPSolver & solver,
			double & min_linear_problem, double & linear_cost_b, double & overall_cost_b, DMatrix & b) const;

	void solve_linear_problem_(const QAPInstance_ & qap_instance, LSAPESolver & solver,
			double & min_linear_problem, double & linear_cost_b, double & overall_cost_b, DMatrix & b) const;

	void solver_to_matrix_(const LSAPSolver & solver, DMatrix & b) const;

	void solver_to_matrix_(const LSAPESolver & solver, DMatrix & b) const;

	void init_linear_cost_matrix_(const QAPInstance_ & qap_instance, DMatrix & linear_cost_matrix) const;

	void init_next_linear_problem_(const QAPInstance_ & qap_instance, const DMatrix & x, const DMatrix & linear_cost_matrix, DMatrix & linear_problem) const;

	bool termination_criterion_met_(const Timer & timer, const double & alpha, const double & min_linear_problem, const std::size_t & current_itr, double lower_bound, double upper_bound) const;

};

}

#endif /* SRC_METHODS_IPFP_HPP_ */
