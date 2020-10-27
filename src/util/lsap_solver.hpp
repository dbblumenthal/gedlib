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
 * @file  lsap_solver.hpp
 * @brief ged::LSAPSolver class declaration.
 */

#ifndef SRC_UTIL_LSAP_SOLVER_HPP_
#define SRC_UTIL_LSAP_SOLVER_HPP_

#include "../env/common_types.hpp"
#include "../env/matrix.hpp"

namespace ged {

/*!
 * @brief This class solves LSAP instances by calling the library lsape available at https://bougleux.users.greyc.fr/lsape/.
 */
class LSAPSolver {

public:

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
	 * @brief Constructs solver for LSAP problem instance.
	 * @param[in] cost_matrix Pointer to the LSAP problem instance that should be solved.
	 */
	LSAPSolver(const DMatrix * cost_matrix);

	/*!
	 * @brief Constructs empty solver.
	 */
	LSAPSolver();

	/*!
	 * @brief Solves the LSAP problem instance.
	 * @param[in] num_solutions The maximal number of solutions that should be computed.
	 */
	void solve(int num_solutions = 1);

	/*!
	 * @brief Sets the LSAP problem instance.
	 * @param[in] cost_matrix Pointer to the LSAP problem instance that should be solved.
	 */
	void set_problem(const DMatrix * cost_matrix);

	/*!
	 * @brief Makes the solver use a greedy method.
	 * @param[in] greedy_method The greedy method that should be sued.
	 */
	void set_greedy_method(const GreedyMethod & greedy_method);

	/*!
	 * @brief Specifies the method used for enumerating multiple solutions.
	 * @param[in] enumeration_method The enumeration method that should be used.
	 */
	void set_enumeration_method(const EnumerationMethod & enumeration_method);

	/*!
	 * @brief Makes the solver use the Hungarian algorithm for optimal solving.
	 */
	void set_hungarian_algorithm();

	/*!
	 * @brief Clears a previously computed solution.
	 */
	void clear_solution();

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
	 * @brief Returns the LSAP problem instance.
	 * @returns Constant pointer to the LSAP problem instance.
	 */
	const DMatrix * cost_matrix() const;

	/*!
	 * @brief Returns the number of rows of the LSAP problem instance.
	 * @returns Number of rows of the LSAP problem instance.
	 */
	std::size_t num_rows() const;

	/*!
	 * @brief Returns the number of columns of the LSAP problem instance.
	 * @returns Number of columns of the LSAP problem instance.
	 */
	std::size_t num_cols() const;

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
	 * @brief Returns the number of solutions.
	 * @returns Actual number of solutions computed by solve(). Might be smaller than @p num_solutions.
	 */
	std::size_t num_solutions() const;

private:

	const DMatrix * cost_matrix_;

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

#include "lsap_solver.ipp"

#endif /* SRC_UTIL_LSAP_SOLVER_HPP_ */
