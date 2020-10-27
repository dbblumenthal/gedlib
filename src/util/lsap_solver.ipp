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
 * @file  lsap_solver.ipp
 * @brief ged::LSAPSolver class definition.
 */

#ifndef SRC_UTIL_LSAP_SOLVER_IPP_
#define SRC_UTIL_LSAP_SOLVER_IPP_

namespace ged {

LSAPSolver ::
LSAPSolver(const DMatrix * cost_matrix) :
cost_matrix_{cost_matrix},
greedy_method_{BASIC},
enumeration_method_{DISSIMILAR},
solve_optimally_{true},
minimal_cost_{0.0},
row_to_col_assignments_(),
col_to_row_assignments_(),
dual_var_rows_(num_rows()),
dual_var_cols_(num_cols()) {}

LSAPSolver ::
LSAPSolver() :
cost_matrix_{nullptr},
greedy_method_{BASIC},
enumeration_method_{DISSIMILAR},
solve_optimally_{true},
minimal_cost_{0.0},
row_to_col_assignments_(),
col_to_row_assignments_(),
dual_var_rows_(),
dual_var_cols_() {}

void
LSAPSolver::
set_problem(const DMatrix * cost_matrix) {
	cost_matrix_ = cost_matrix;
	clear_solution();
}

void
LSAPSolver::
set_greedy_method(const GreedyMethod & greedy_method) {
	solve_optimally_ = false;
	greedy_method_ = greedy_method;
}

void
LSAPSolver::
set_enumeration_method(const EnumerationMethod & enumeration_method) {
    enumeration_method_ = enumeration_method;
}

void
LSAPSolver::
set_hungarian_algorithm() {
	solve_optimally_ = true;
}

void
LSAPSolver ::
clear_solution() {
	minimal_cost_ = 0.0;
	row_to_col_assignments_.clear();
	col_to_row_assignments_.clear();
	dual_var_rows_ = std::vector<double>(num_rows());
	dual_var_cols_ = std::vector<double>(num_cols());
}

void
LSAPSolver ::
solve(int num_solutions) {
	clear_solution();
	if (not solve_optimally_) {
	    std::cout << "WARNING: greedy solvers not implemented yet. Using optimal solvers.";
	}
    DMatrix copied_costs(*cost_matrix_);
    liblsap::LSAP<double, int> lsap(static_cast<int>(num_rows()), static_cast<int>(num_cols()), copied_costs.data());
    lsap.solve();
    if (num_solutions > 1) {
        if (enumeration_method_ == BASELINE) {
            lsap.enumerate(static_cast<int>(num_solutions), liblsap::ENUM_EDG_SELECT_BALANCED);
        }
        else {
            lsap.enumerateDissimilar(static_cast<int>(num_solutions), liblsap::ENUM_DISS_MXW_DFS,
                                     liblsap::ENUM_EDG_SELECT_RAND, liblsap::ENUM_CYCLE_SELECT_RAND);
        }
    }
    const std::vector<liblsap::Matching<int>> & matchings{lsap.primals()};
    for (const liblsap::Matching<int> & matching : matchings) {
        row_to_col_assignments_.push_back(std::vector<std::size_t>(num_rows()));
        for (int row{0}; row < static_cast<int>(num_rows()); row++) {
            row_to_col_assignments_.back().at(row) = matching(0, row);
        }
        col_to_row_assignments_.push_back(std::vector<std::size_t>(num_cols()));
        for (int col{0}; col < static_cast<int>(num_cols()); col++) {
            col_to_row_assignments_.back().at(col) = matching(1, col);
        }
    }
    for (int row{0}; row < static_cast<int>(num_rows()); row++) {
        dual_var_rows_.emplace_back(lsap.dual(0, row));
    }
    for (int col{0}; col < static_cast<int>(num_cols()); col++) {
        dual_var_cols_.emplace_back(lsap.dual(1, col));
    }
    compute_cost_from_dual_vars_();


	/* Old implementation.
	if (solve_optimally_) {
		lsape::hungarianLSAP(cost_matrix_->data(), num_rows(), num_cols(), row_to_col_assignments_.at(0).data(), dual_var_rows_.data(), dual_var_cols_.data(), col_to_row_assignments_.at(0).data());
		compute_cost_from_dual_vars_();
		if (num_solutions > 1) {
			std::list<std::size_t *> further_row_to_col_assignments;
			lsape::lsapSolutions(cost_matrix_->data(), num_rows(), num_cols(), num_solutions, row_to_col_assignments_.at(0).data(), dual_var_rows_.data(), dual_var_cols_.data(),
					further_row_to_col_assignments);
			row_to_col_assignments_.clear();
			col_to_row_assignments_.clear();
			for (auto row_to_col_assignment = further_row_to_col_assignments.begin(); row_to_col_assignment != further_row_to_col_assignments.end(); row_to_col_assignment++) {
				row_to_col_assignments_.push_back(std::vector<std::size_t>(*row_to_col_assignment, *row_to_col_assignment + num_rows()));
				col_to_row_assignments_.push_back(construct_col_to_row_asignment_(row_to_col_assignments_.back()));
				delete *row_to_col_assignment;
			}
		}
	}
	else {
		minimal_cost_ = lsape::greedyLSAP(cost_matrix_->data(), num_rows(), num_cols(), row_to_col_assignments_.at(0).data(), col_to_row_assignments_.at(0).data(), greedy_method_);
	}
	*/
}

double
LSAPSolver ::
minimal_cost() const {
	return minimal_cost_;
}

std::size_t
LSAPSolver ::
get_assigned_col(std::size_t row, std::size_t solution_id) const {
	return row_to_col_assignments_.at(solution_id).at(row);
}

std::size_t
LSAPSolver ::
get_assigned_row(std::size_t col, std::size_t solution_id) const {
	return col_to_row_assignments_.at(solution_id).at(col);
}

double
LSAPSolver ::
get_slack(std::size_t row, std::size_t col) const {
	return cost_matrix_->operator()(row, col) - (dual_var_rows_.at(row) + dual_var_cols_.at(col));
}

double
LSAPSolver ::
get_dual_var_row(std::size_t row) const {
	return dual_var_rows_.at(row);
}

double
LSAPSolver ::
get_dual_var_col(std::size_t col) const {
	return dual_var_cols_.at(col);
}

const std::vector<std::size_t> &
LSAPSolver ::
row_to_col_assignment(std::size_t solution_id) const {
	return row_to_col_assignments_.at(solution_id);
}

const std::vector<std::size_t> &
LSAPSolver ::
col_to_row_assignment(std::size_t solution_id) const {
	return col_to_row_assignments_.at(solution_id);
}

const DMatrix *
LSAPSolver ::
cost_matrix() const {
	return cost_matrix_;
}

std::size_t
LSAPSolver ::
num_rows() const {
	return cost_matrix_->num_rows();
}

std::size_t
LSAPSolver ::
num_cols() const {
	return cost_matrix_->num_cols();
}

std::size_t
LSAPSolver ::
num_solutions() const {
	return row_to_col_assignments_.size();
}

std::vector<std::size_t>
LSAPSolver ::
construct_col_to_row_asignment_(const std::vector<std::size_t> row_to_col_assignment) const {
	std::vector<std::size_t> col_to_row_assignment(num_cols(), num_rows());
	for (std::size_t row{}; row < num_rows(); row++) {
		if (row_to_col_assignment.at(row) < num_cols()) {
			col_to_row_assignment[row_to_col_assignment.at(row)] = row;
		}
	}
	return col_to_row_assignment;
}

void
LSAPSolver ::
compute_cost_from_assignments_() {
	minimal_cost_ = 0.0;
	for (std::size_t row{}; row < num_rows(); row++) {
		minimal_cost_ += cost_matrix_->operator()(row, row_to_col_assignments_.at(0).at(row));
	}
	for (std::size_t col{}; col < num_cols(); col++) {
		if (col_to_row_assignments_.at(0).at(col) == num_rows()) {
			minimal_cost_ += cost_matrix_->operator()(num_rows(), col);
		}
	}
}

void
LSAPSolver ::
compute_cost_from_dual_vars_() {
	minimal_cost_ = 0.0;
	for (auto dual_var = dual_var_rows_.begin(); dual_var != dual_var_rows_.end(); dual_var++) {
		minimal_cost_ += *dual_var;
	}
	for (auto dual_var = dual_var_cols_.begin(); dual_var != dual_var_cols_.end(); dual_var++)  {
		minimal_cost_ += *dual_var;
	}
}

}

#endif /* SRC_UTIL_LSAP_SOLVER_IPP_ */
