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
 * @file  lsape_solver.hpp
 * @brief ged::LSAPESolver class declaration.
 */

#ifndef LSAPE_SOLVER_HPP
#define LSAPE_SOLVER_HPP 1

#include "../env/common_types.hpp"
#include "../env/matrix.hpp"

namespace ged {

/*!
 * @brief This class solves LSAPE instances by calling the library lsape available at https://bougleux.users.greyc.fr/lsape/.
 */
class LSAPESolver {

public:

	/*!
	 * @brief Selects a model for solving LSAPE with the Hungarian algorithm.
	 * @details The different models are described in:
	 * - S. Bougleux, B. Ga&uuml;z&egrave;re, D. B. Blumenthal, and L. Brun:
	 *   &ldquo;Fast linear sum assignment with error-correction and no cost constraints&rdquo;,
	 *   https://doi.org/10.1016/j.patrec.2018.03.032
	 */
	enum Model {
		ECBP=0, //!< Adaption of Hungarian Algorithm to LSAPE.
		FLWC=1, //!< Reduction to compact LSAP instance without cost constraints.
		EBP=2,  //!< Reduction to extended LSAP instance without cost constraints.
		FLCC=3, //!< Reduction to compact LSAP instance for quasimetric costs.
		FBP=4,  //!< Reduction to compact LSAP instance for quasimetric costs.
		FBP0=5, //!< Reduction to compact LSAP instance for quasimetric costs.
		SFBP=6  //!< Reduction to compact LSAP instance for quasimetric costs.
	};

	/*!
	 * @brief Selects a greedy method.
	 * @details The different greedy methods are described in:
	 * - K. Riesen, M. Ferrer, A. Fischer, and H. Bunke:
	 *   &ldquo;Approximation of graph edit distance in quadratic time&rdquo;,
	 *   https://doi.org/10.1007/978-3-319-18224-7_1
	 */
	enum GreedyMethod {
		BASIC=0,         //!< See https://doi.org/10.1007/978-3-319-18224-7_1.
		REFINED=1,		 //!< See https://doi.org/10.1007/978-3-319-18224-7_1.
		LOSS=2,			 //!< See https://doi.org/10.1007/978-3-319-18224-7_1.
		BASIC_SORT=3,	 //!< See https://doi.org/10.1007/978-3-319-18224-7_1.
		INT_BASIC_SORT=4 //!< See https://doi.org/10.1007/978-3-319-18224-7_1.
	};

	enum EnumerationMethod {
	    BASELINE=0,
	    DISSIMILAR=1
	};

	/*!
	 * @brief Constructs solver for LSAPE problem instance.
	 * @param[in] cost_matrix Pointer to the LSAPE problem instance that should be solved.
	 */
	LSAPESolver(const DMatrix * cost_matrix);

	/*!
	 * @brief Constructs empty solver.
	 */
	LSAPESolver();

	/*!
	 * @brief Sets the LSAPE problem instance.
	 * @param[in] cost_matrix Pointer to the LSAPE problem instance that should be solved.
	 */
	void set_problem(const DMatrix * cost_matrix);

	/*!
	 * @brief Clears a previously computed solution.
	 */
	void clear_solution();

	/*!
	 * @brief Makes the solver use a specific model for optimal solving.
	 * @param[in] model The model that should be used.
	 */
	void set_model(const Model & model);

	/*!
	 * @brief Makes the solver use a greedy method.
	 * @param[in] greedy_method The greedy method that should be used.
	 */
	void set_greedy_method(const GreedyMethod & greedy_method);

    /*!
     * @brief Specifies the method used for enumerating multiple solutions.
     * @param[in] enumeration_method The enumeration method that should be used.
     */
    void set_enumeration_method(const EnumerationMethod & enumeration_method);

	/*!
	 * @brief Solves the LSAPE problem instance.
	 * @param[in] num_solutions The maximal number of solutions that should be computed.
	 */
	void solve(int num_solutions = 1);

	/*!
	 * @brief Returns the cost of the computed solutions.
	 * @return Cost of computed solutions.
	 */
	double minimal_cost() const;

	/*!
	 * @brief Returns the assigned column.
	 * @param[in] row Row whose assigned column should be returned.
	 * @param[in] solution_id ID of the solution where the assignment should be looked up.
	 * @returns Column to which @p row is assigned to in solution with ID @p solution_id or ged::undefined() if @p row is not assigned to any column.
	 */
	std::size_t get_assigned_col(std::size_t row, std::size_t solution_id = 0) const;

	/*!
	 * @brief Returns the assigned row.
	 * @param[in] col Column whose assigned row should be returned.
	 * @param[in] solution_id ID of the solution where the assignment should be looked up.
	 * @returns Row to which @p col is assigned to in solution with ID @p solution_id or ged::undefined() if @p col is not assigned to any row.
	 */
	std::size_t get_assigned_row(std::size_t col, std::size_t solution_id = 0) const;

	/*!
	 * @brief Returns the slack of a cell.
	 * @param[in] row Row of the cell.
	 * @param[in] col Column of the cell.
	 * @returns Slack of the cell (@p row, @p col).
	 */
	double get_slack(std::size_t row, std::size_t col) const;

	/*!
	 * @brief Returns the dual variable of a row.
	 * @param[in] row The row.
	 * @returns Dual variable of the row @p row.
	 */
	double get_dual_var_row(std::size_t row) const;

	/*!
	 * @brief Returns the dual variable of a column.
	 * @param[in] col The column.
	 * @returns Dual variable of the column @p col.
	 */
	double get_dual_var_col(std::size_t col) const;

	/*!
	 * @brief Returns the LSAPE problem instance.
	 * @returns Constant pointer to the LSAPE problem instance.
	 */
	const DMatrix * cost_matrix() const;

	/*!
	 * @brief Returns the assignment from rows to columns.
	 * @param[in] solution_id The ID of the solutions whose assignment should be returned.
	 * @returns Constant reference of assignment of rows to columns in the solution with ID @p solution_id.
	 */
	const std::vector<std::size_t> & row_to_col_assignment(std::size_t solution_id = 0) const;

	/*!
	 * @brief Returns the assignment from columns to rows.
	 * @param[in] solution_id The ID of the solutions whose assignment should be returned.
	 * @returns Constant reference of assignment of columns to rows in the solution with ID @p solution_id.
	 */
	const std::vector<std::size_t> & col_to_row_assignment(std::size_t solution_id = 0) const;

	/*!
	 * @brief Returns the number of real rows of the LSAPE problem instance, i.e, the number of rows of the internal cost matrix minus 1.
	 * @returns Number of real rows of the LSAPE problem instance.
	 */
	std::size_t num_rows() const;

	/*!
	 * @brief Returns the number of real columns of the LSAPE problem instance, i.e, the number of columns of the internal cost matrix minus 1.
	 * @returns Number of real columns of the LSAPE problem instance.
	 */
	std::size_t num_cols() const;

	/*!
	 * @brief Returns the total number of rows of the LSAPE problem instance.
	 * @returns Total number of rows of the LSAPE problem instance.
	 */
	std::size_t total_num_rows() const;

	/*!
	 * @brief Returns the total number of columns of the LSAPE problem instance.
	 * @returns Total number of columns of the LSAPE problem instance.
	 */
	std::size_t total_num_cols() const;

	/*!
	 * @brief Returns the number of solutions.
	 * @returns Actual number of solutions computed by solve(). Might be smaller than @p num_solutions.
	 */
	std::size_t num_solutions() const;


private:

	const DMatrix * cost_matrix_;

	Model model_;

	GreedyMethod greedy_method_;

	EnumerationMethod enumeration_method_;

	bool solve_optimally_;

	double minimal_cost_;

	std::vector<std::vector<std::size_t>> row_to_col_assignments_;

	std::vector<std::vector<std::size_t>> col_to_row_assignments_;

	std::vector<double> dual_var_rows_;

	std::vector<double> dual_var_cols_;

	// Helper functions.
	std::vector<std::size_t> construct_col_to_row_asignment_(const std::vector<std::size_t> row_to_col_assignment) const;

	void compute_cost_from_assignments_();

	void compute_cost_from_dual_vars_();

};

}

#include "lsape_solver.ipp"

#endif
