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
 * \file  lsape_solver.ipp
 * \brief ged::LSAPESolver class definition.
 */

#ifndef LSAPE_SOLVER_IPP
#define LSAPE_SOLVER_IPP 1

namespace ged {

LSAPESolver ::
LSAPESolver(const DMatrix * cost_matrix) :
cost_matrix_{cost_matrix},
model_{ECBP},
greedy_method_{BASIC},
enumeration_method_{DISSIMILAR},
solve_optimally_{true},
minimal_cost_{0.0},
row_to_col_assignments_(),
col_to_row_assignments_(),
dual_var_rows_(total_num_rows()),
dual_var_cols_(total_num_cols()) {}

LSAPESolver ::
LSAPESolver() :
cost_matrix_{nullptr},
model_{ECBP},
greedy_method_{BASIC},
enumeration_method_{DISSIMILAR},
solve_optimally_{true},
minimal_cost_{0.0},
row_to_col_assignments_(),
col_to_row_assignments_() {}

void
LSAPESolver::
set_problem(const DMatrix * cost_matrix) {
	cost_matrix_ = cost_matrix;
}

void
LSAPESolver ::
clear_solution() {
	minimal_cost_ = 0.0;
	row_to_col_assignments_.clear();
	col_to_row_assignments_.clear();
	dual_var_rows_ = std::vector<double>(total_num_rows());
	dual_var_cols_ = std::vector<double>(total_num_cols());
}

void
LSAPESolver::
set_model(const Model & model) {
	solve_optimally_ = true;
	model_ = model;
}

void
LSAPESolver::
set_greedy_method(const GreedyMethod & greedy_method) {
	solve_optimally_ = false;
	greedy_method_ = greedy_method;
}

void
LSAPESolver::
set_enumeration_method(const EnumerationMethod &enumeration_method) {
    enumeration_method_ = enumeration_method;
}

void
LSAPESolver ::
solve(int num_solutions) {
	clear_solution();

    if (not solve_optimally_) {
        std::cout << "WARNING: greedy solvers not implemented yet. Using optimal solvers.";
    }
    DMatrix copied_costs(*cost_matrix_);
    liblsap::LSAPE<double, int> lsape(static_cast<int>(num_rows()), static_cast<int>(num_cols()), copied_costs.data());
    lsape.arithScaleInstance();
   // std::string instance_file("ged_instance_" + std::to_string(num_rows()) + "_" + std::to_string(num_cols()) + ".txt"); // todo rm
   // std::cout << "saving instance " << instance_file << std::endl; // todo rm
   // lsape.saveInstance(instance_file.c_str()); // todo rm
    lsape.solve();
    if (num_solutions > 1) {
        if (enumeration_method_ == BASELINE) {
            lsape.enumerate(static_cast<int>(num_solutions), liblsap::ENUM_EDG_SELECT_BALANCED);
        }
        else {
            lsape.enumerateDissimilar(static_cast<int>(num_solutions), liblsap::ENUM_DISS_MXW_DFS_DUMMY1_FREE,liblsap::ENUM_EDG_SELECT_RAND, liblsap::ENUM_CYCLE_SELECT_RAND);
        }
    }
    const std::vector<liblsap::Matching<int>> & matchings{lsape.primals()};
    for (const liblsap::Matching<int> & matching : matchings) {
        row_to_col_assignments_.push_back(std::vector<std::size_t>(num_rows()));
        for (int row{0}; row < static_cast<int>(num_rows()); row++) {
            row_to_col_assignments_.back()[row] = matching(0, row);
        }
        col_to_row_assignments_.push_back(std::vector<std::size_t>(num_cols()));
        for (int col{0}; col < static_cast<int>(num_cols()); col++) {
            col_to_row_assignments_.back()[col] = matching(1, col);
        }
    }
    for (int row{0}; row < static_cast<int>(total_num_rows()); row++) {
        dual_var_rows_.emplace_back(lsape.dual(0, row));
    }
    for (int col{0}; col < static_cast<int>(total_num_cols()); col++) {
        dual_var_cols_.emplace_back(lsape.dual(1, col));
    }
    compute_cost_from_dual_vars_();
	/* Old implementation.
	if (solve_optimally_) {
		lsape::lsapeSolver(cost_matrix_->data(), total_num_rows(), total_num_cols(), row_to_col_assignments_.at(0).data(), col_to_row_assignments_.at(0).data(),
				dual_var_rows_.data(), dual_var_cols_.data(), static_cast<lsape::LSAPE_MODEL>(model_));
		compute_cost_from_assignments_();
		if (num_solutions > 1) {
			std::list<std::size_t *> further_row_to_col_assignments;
			lsape::lsapeSolutionsFromOne(cost_matrix_->data(), total_num_rows(),total_num_cols(), num_solutions, row_to_col_assignments_.at(0).data(), col_to_row_assignments_.at(0).data(),
					dual_var_rows_.data(), dual_var_cols_.data(), further_row_to_col_assignments);
			row_to_col_assignments_.clear();
			col_to_row_assignments_.clear();
			for (auto row_to_col_assignment : further_row_to_col_assignments) {
				row_to_col_assignments_.push_back(std::vector<std::size_t>(row_to_col_assignment, row_to_col_assignment + num_rows()));
				col_to_row_assignments_.push_back(construct_col_to_row_asignment_(row_to_col_assignments_.back()));
			}
			for (auto && row_to_col_assignment : further_row_to_col_assignments) {
				delete row_to_col_assignment;
			}
		}
	}
	else {
		minimal_cost_ = lsape::lsapeGreedy(cost_matrix_->data(), num_rows(), num_cols(), row_to_col_assignments_.at(0).data(), col_to_row_assignments_.at(0).data(),
				static_cast<lsape::GREEDY_METHOD>(greedy_method_));
	}
	*/
}

double
LSAPESolver ::
minimal_cost() const {
	return minimal_cost_;
}

std::size_t
LSAPESolver ::
get_assigned_col(LabelID row, std::size_t solution_id) const {
	return row_to_col_assignments_.at(solution_id).at(row);
}

std::size_t
LSAPESolver ::
get_assigned_row(LabelID col, std::size_t solution_id) const {
	return col_to_row_assignments_.at(solution_id).at(col);
}

const std::vector<std::size_t> &
LSAPESolver ::
row_to_col_assignment(std::size_t solution_id) const {
	return row_to_col_assignments_.at(solution_id);
}

const std::vector<std::size_t> &
LSAPESolver ::
col_to_row_assignment(std::size_t solution_id) const {
	return col_to_row_assignments_.at(solution_id);
}

double
LSAPESolver ::
get_slack(std::size_t row, std::size_t col) const {
	return cost_matrix_->operator()(row, col) - (dual_var_rows_.at(row) + dual_var_cols_.at(col));
}

double
LSAPESolver ::
get_dual_var_row(std::size_t row) const {
	return dual_var_rows_.at(row);
}

double
LSAPESolver ::
get_dual_var_col(std::size_t col) const {
	return dual_var_cols_.at(col);
}

const DMatrix *
LSAPESolver ::
cost_matrix() const {
	return cost_matrix_;
}

std::size_t
LSAPESolver ::
num_rows() const {
	return cost_matrix_->num_rows() - 1;
}

std::size_t
LSAPESolver ::
num_cols() const {
	return cost_matrix_->num_cols() - 1;
}

std::size_t
LSAPESolver ::
total_num_rows() const {
	return cost_matrix_->num_rows();
}

std::size_t
LSAPESolver ::
total_num_cols() const {
	return cost_matrix_->num_cols();
}

std::size_t
LSAPESolver ::
num_solutions() const {
	return row_to_col_assignments_.size();
}

std::vector<std::size_t>
LSAPESolver ::
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
LSAPESolver ::
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
LSAPESolver ::
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


#endif
