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

#include "catch.hpp"
#include <iomanip>

#include "../../../src/util/lsap_solver.hpp"
#include "../../../src/util/lsape_solver.hpp"

std::ostream& bold_on(std::ostream& os)
{
	return os << "\033[1m";
}

std::ostream& bold_off(std::ostream& os)
{
	return os << "\033[0m";
}

void print_solution(const ged::LSAPESolver & solver, const std::string heading = "") {
	std::cout << "\t" << heading << std::endl << "\tdual variables and slack graph" << std::endl << "\t   " << bold_on;
	for (std::size_t col{}; col < solver.total_num_cols(); col++) {
		std::cout << std::setw(2) << std::right << solver.get_dual_var_col(col) << " ";
	}
	std::cout << bold_off;
	for (std::size_t row{}; row < solver.total_num_rows(); row++) {
		std::cout << std::endl << "\t" << bold_on << std::setw(2) << std::right << solver.get_dual_var_row(row) << bold_off << " ";
		for (std::size_t col{}; col < solver.total_num_cols(); col++) {
			std::cout << std::setw(2) << std::right << solver.get_slack(row, col) << " ";
		}
	}
	std::cout << std::endl;
	for (std::size_t solution_id{0}; solution_id < solver.num_solutions(); solution_id++) {
		std::vector<std::size_t> substituted_rows;
		std::vector<std::size_t> substituted_cols;
		std::vector<std::size_t> deleted_rows;
		std::vector<std::size_t> inserted_cols;
		for (std::size_t row{}; row < solver.num_rows(); row++) {
			if (solver.get_assigned_col(row, solution_id) < solver.num_cols()) {
				substituted_rows.push_back(row);
			}
			else {
				deleted_rows.push_back(row);
			}
		}
		for (std::size_t col{}; col < solver.num_cols(); col++) {
			if (solver.get_assigned_row(col, solution_id) < solver.num_rows()) {
				substituted_cols.push_back(col);
			}
			else {
				inserted_cols.push_back(col);
			}
		}
		std::cout << "\tsolution_id = " << solution_id << std::endl;
		std::cout << "\t\tsubstitutions (row->col)" << std::endl << "\t\t{";
		for (auto row : substituted_rows) {
			std::cout << " " << row << "->" << solver.get_assigned_col(row, solution_id);
		}
		std::cout << " }" << std::endl << "\t\tdeleted rows" << std::endl << "\t\t{";
		for (auto row : deleted_rows) {
			std::cout << " " << row;
		}
		std::cout << " }" << std::endl << "\t\tinserted columns" << std::endl << "\t\t{";
		for (auto col : inserted_cols) {
			std::cout << " " << col;
		}
		std::cout << " }" << std::endl;
	}
}

void print_solution(const ged::LSAPSolver & solver, const std::string heading = "") {
	std::cout << "\t" << heading << std::endl << "\tdual variables and slack graph" << std::endl << "\t   " << bold_on;
	for (std::size_t col{}; col < solver.num_cols(); col++) {
		std::cout << std::setw(2) << std::right << solver.get_dual_var_col(col) << " ";
	}
	std::cout << bold_off;
	for (std::size_t row{}; row < solver.num_rows(); row++) {
		std::cout << std::endl << "\t" << bold_on << std::setw(2) << std::right << solver.get_dual_var_row(row) << bold_off << " ";
		for (std::size_t col{}; col < solver.num_cols(); col++) {
			std::cout << std::setw(2) << std::right << solver.get_slack(row, col) << " ";
		}
	}
	std::cout << std::endl;
	for (std::size_t solution_id{0}; solution_id < solver.num_solutions(); solution_id++) {
		std::vector<std::size_t> substituted_rows;
		std::vector<std::size_t> substituted_cols;
		std::vector<std::size_t> deleted_rows;
		std::vector<std::size_t> inserted_cols;
		for (std::size_t row{}; row < solver.num_rows(); row++) {
			if (solver.get_assigned_col(row, solution_id) < solver.num_cols()) {
				substituted_rows.push_back(row);
			}
			else {
				deleted_rows.push_back(row);
			}
		}
		for (std::size_t col{}; col < solver.num_cols(); col++) {
			if (solver.get_assigned_row(col, solution_id) < solver.num_rows()) {
				substituted_cols.push_back(col);
			}
			else {
				inserted_cols.push_back(col);
			}
		}
		std::cout << "\tsolution_id = " << solution_id << std::endl;
		std::cout << "\t\tsubstitutions (row->col)" << std::endl << "\t\t{";
		for (auto row : substituted_rows) {
			std::cout << " " << row << "->" << solver.get_assigned_col(row, solution_id);
		}
		std::cout << " }" << std::endl << "\t\tdeleted rows" << std::endl << "\t\t{";
		for (auto row : deleted_rows) {
			std::cout << " " << row;
		}
		std::cout << " }" << std::endl << "\t\tinserted collumns" << std::endl << "\t\t{";
		for (auto col : inserted_cols) {
			std::cout << " " << col;
		}
		std::cout << " }" << std::endl;
	}
}

TEST_CASE("LSAPESolver works") {

	SECTION("multiple solutions 1") {
		ged::DMatrix ecm_mult_1(3,3);
		ecm_mult_1.matrix() <<
				9, 2, 3,
				2, 9, 1,
				1, 3, 0;
		std::cout << std::endl << "===multiple solutions LSAPE 1===" << std::endl;
		std::cout << "LSAPE instance" << std::endl << ecm_mult_1.matrix() << std::endl;
		ged::LSAPESolver solver_mult_sol;
		solver_mult_sol.set_problem(&ecm_mult_1);
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FLWC);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FLWC multiple solutions 1===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FLCC);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FLCC multiple solutions 1===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::ECBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===ECBP multiple solutions 1===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::EBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===EBP multiple solutions 1===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FBP multiple solutions 1===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::SFBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===SFBP multiple solutions 1===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FBP0);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FBP0 multiple solutions 1===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
	}

	SECTION("multiple solutions 2") {
		ged::DMatrix ecm_mult_2(2,2);
		ecm_mult_2.matrix() <<
				2, 1,
				1, 0;
		std::cout << std::endl << "===multiple solutions LSAPE 2===" << std::endl;
		std::cout << "LSAPE instance" << std::endl << ecm_mult_2.matrix() << std::endl;
		ged::LSAPESolver solver_mult_sol;
		solver_mult_sol.set_problem(&ecm_mult_2);
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FLWC);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FLWC multiple solutions 2===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(2.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FLCC);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FLCC multiple solutions 2===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(2.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::ECBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===ECBP multiple solutions 2===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(2.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::EBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===EBP multiple solutions 2===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(2.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FBP multiple solutions 2===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(2.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::SFBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===SFBP multiple solutions 2===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(2.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FBP0);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FBP0 multiple solutions 2===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(2.0));
	}

	SECTION("multiple solutions 3") {
		ged::DMatrix ecm_mult_3(3,3);
		ecm_mult_3.matrix() <<
				9, 3, 1,
				2, 9, 1,
				1, 1, 0;
		std::cout << std::endl << "===multiple solutions LSAPE 3===" << std::endl;
		std::cout << "LSAPE instance" << std::endl << ecm_mult_3.matrix() << std::endl;
		ged::LSAPESolver solver_mult_sol;
		solver_mult_sol.set_problem(&ecm_mult_3);
		solver_mult_sol.set_model(ged::LSAPESolver::Model::FLWC);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===FLWC multiple solutions 3===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::ECBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===ECBP multiple solutions 3===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
		solver_mult_sol.set_model(ged::LSAPESolver::Model::EBP);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===EBP multiple solutions 3===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(4.0));
	}

	ged::DMatrix ecm1(5,6);

	ecm1.matrix() <<
			7, 11, 9, 8, 9, 10,
			2,  8, 8, 5, 7,  3,
			1,  7, 6, 6, 9,  5,
			3,  7, 6, 2, 2,  3,
			4,  2, 2, 7, 8,  0;


	SECTION("matrix setup correctly") {
		CHECK(ecm1(1,2) == 8);
		CHECK(ecm1(3,0) == 3);
		CHECK(ecm1(0,0) == 7);
		CHECK(ecm1(4,5) == 0);
	}
	std::size_t num_rows{4};
	std::size_t num_cols{5};
	ged::LSAPESolver solver;
	solver.set_problem(&ecm1);
	solver.set_model(ged::LSAPESolver::Model::ECBP);
	solver.solve(10);
	double true_min{solver.minimal_cost()};

	SECTION("HNG works as expected") {
		//print_solution(solver, "===ECBP===");
		for (std::size_t row{}; row < num_rows; row++)
			if (solver.get_assigned_col(row) < num_cols)
				CHECK(solver.get_assigned_row(solver.get_assigned_col(row)) == row);
		for (std::size_t col{}; col < num_cols; col++)
			if (solver.get_assigned_row(col) < num_rows)
				CHECK(solver.get_assigned_col(solver.get_assigned_row(col)) == col);
	}

	SECTION("EBP works") {
		solver.set_model(ged::LSAPESolver::Model::EBP);
		solver.solve(10);
		//print_solution(solver, "===EBP===");
		for (std::size_t row{}; row < num_rows; row++)
			if (solver.get_assigned_col(row) < num_cols)
				CHECK(solver.get_assigned_row(solver.get_assigned_col(row)) == row);
		for (std::size_t col{}; col < num_cols; col++)
			if (solver.get_assigned_row(col) < num_rows)
				CHECK(solver.get_assigned_col(solver.get_assigned_row(col)) == col);
		CHECK(true_min == Approx(solver.minimal_cost()));
	}
	SECTION("FLWC works") {
		solver.set_model(ged::LSAPESolver::Model::FLWC);
		solver.solve(10);
		//print_solution(solver, "===FLWC===");
		for (std::size_t row{}; row < num_rows; row++)
			if (solver.get_assigned_col(row) < num_cols)
				CHECK(solver.get_assigned_row(solver.get_assigned_col(row)) == row);
		for (std::size_t col{}; col < num_cols; col++)
			if (solver.get_assigned_row(col) < num_rows)
				CHECK(solver.get_assigned_col(solver.get_assigned_row(col)) == col);
		CHECK(true_min == Approx(solver.minimal_cost()));
	}
	SECTION("FBP fails as exppected (trianlge-ineq not satisfied)") {
		solver.set_model(ged::LSAPESolver::Model::FBP);
		solver.solve(10);
		//print_solution(solver, "===FBP===");
		for (std::size_t row{}; row < num_rows; row++)
			if (solver.get_assigned_col(row) < num_cols)
				CHECK(solver.get_assigned_row(solver.get_assigned_col(row)) == row);
		for (std::size_t col{}; col < num_cols; col++)
			if (solver.get_assigned_row(col) < num_rows)
				CHECK(solver.get_assigned_col(solver.get_assigned_row(col)) == col);
		CHECK(true_min <= Approx(solver.minimal_cost()));
	}
	SECTION("SFBP fails as exppected (trianlge-ineq not satisfied)") {
		solver.set_model(ged::LSAPESolver::Model::SFBP);
		solver.solve(10);
		//print_solution(solver, "===SFBP===");
		for (std::size_t row{}; row < num_rows; row++)
			if (solver.get_assigned_col(row) < num_cols)
				CHECK(solver.get_assigned_row(solver.get_assigned_col(row)) == row);
		for (std::size_t col{}; col < num_cols; col++)
			if (solver.get_assigned_row(col) < num_rows)
				CHECK(solver.get_assigned_col(solver.get_assigned_row(col)) == col);
		CHECK(true_min <= Approx(solver.minimal_cost()));
	}
	SECTION("FBP0 fails as exppected (trianlge-ineq not satisfied)") {
		solver.set_model(ged::LSAPESolver::Model::FBP0);
		solver.solve(10);
		//print_solution(solver, "===FBP0===");
		for (std::size_t row{}; row < num_rows; row++)
			if (solver.get_assigned_col(row) < num_cols)
				CHECK(solver.get_assigned_row(solver.get_assigned_col(row)) == row);
		for (std::size_t col{}; col < num_cols; col++)
			if (solver.get_assigned_row(col) < num_rows)
				CHECK(solver.get_assigned_col(solver.get_assigned_row(col)) == col);
		CHECK(true_min <= Approx(solver.minimal_cost()));
	}

}

TEST_CASE("LSAPSolver works") {

	SECTION("multiple solutions") {
		ged::DMatrix ecm(2,2);
		ecm.matrix() <<
				1, 1,
				2, 2;
		std::cout << std::endl << "===multiple solutions LSAP===" << std::endl;
		std::cout << "LSAP instance" << std::endl << ecm.matrix() << std::endl;
		ged::LSAPSolver solver_mult_sol;
		solver_mult_sol.set_problem(&ecm);
		solver_mult_sol.solve(10);
		print_solution(solver_mult_sol, "===multiple solutions===");
		CHECK(solver_mult_sol.minimal_cost() == Approx(3.0));
	}

	SECTION("Squared matrix") {
		ged::DMatrix matrix_squared(3,3,2.0);
		for (std::size_t row{}; row < 3; ++row) {
			matrix_squared(row, row) = 1.0;
		}

		ged::LSAPSolver solver_squared;
		solver_squared.set_problem(&matrix_squared);
		solver_squared.solve();
		CHECK(solver_squared.minimal_cost() == Approx(3.0));
	}

	SECTION("Matrix with more rows") {

		ged::DMatrix matrix_more_rows(3,2,1.0);
		for (std::size_t row{}; row < 2; ++row) {
			matrix_more_rows(row, row) = 0.0;
		}


		ged::LSAPSolver solver_more_rows;
		solver_more_rows.set_problem(&matrix_more_rows);
		solver_more_rows.solve();
		CHECK(solver_more_rows.minimal_cost() == Approx(0.0));
	}

	SECTION("Matrix with more columns") {

		ged::DMatrix matrix_more_cols_1(3,2,0.0);
		for (std::size_t row = 0; row < 2; ++row) {
			matrix_more_cols_1(row, row) =-1.0;
		}
		ged::LSAPSolver solver_more_cols;
		solver_more_cols.set_problem(&matrix_more_cols_1);
		solver_more_cols.solve();
		CHECK(solver_more_cols.minimal_cost() == Approx(-2.0));

		ged::DMatrix matrix_more_cols_2(2,4,0.0);
		for (std::size_t row = 0; row < 2; ++row) {
			matrix_more_cols_2(row, row) = -1.0;
		}
		solver_more_cols.set_problem(&matrix_more_cols_2);
		solver_more_cols.solve();
		CHECK(solver_more_cols.minimal_cost() == Approx(-2.0));
	}
}




