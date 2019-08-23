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
 * @file progress_bar.ipp
 * @brief ged::ProgressBar class definition.
 */

#ifndef SRC_ENV_PROGRESS_BAR_IPP_
#define SRC_ENV_PROGRESS_BAR_IPP_

namespace ged {

ProgressBar ::
ProgressBar(std::size_t num_tasks) :
num_solved_tasks_{0},
num_tasks_{num_tasks},
start_time_{std::chrono::high_resolution_clock::now()} {}

void
ProgressBar ::
increment() {
	num_solved_tasks_++;
}

void
ProgressBar ::
reset() {
	num_solved_tasks_ = 0;
}

std::ostream &
operator<<(std::ostream & os, const ProgressBar & progress_bar) {

	std::streamsize precision{os.precision()};
	os.precision(2);
	os.setf(std::ios::fixed, std::ios::floatfield);
	double progress_in_percent{100.0 * static_cast<double>(progress_bar.num_solved_tasks_) / static_cast<double>(progress_bar.num_tasks_)};
	os << "[";
	for (std::size_t i{1}; i <= 10; i++) {
		if (static_cast<double>(i) <= progress_in_percent / 10.0) {
			os << "=";
		}
		else if (static_cast<double>(i - 1) < progress_in_percent / 10.0) {
			os << ">";
		}
		else {
			os << " ";
		}
	}
	os << "] ";
	if (progress_in_percent < 100) {
		os << " ";
	}
	if (progress_in_percent < 10) {
		os << " ";
	}
	os << progress_in_percent << " %";
	if (progress_in_percent > 0) {
		std::size_t max_num_digits{7};
		std::size_t num_digits{4};
		std::chrono::duration<double> runtime_so_far = std::chrono::high_resolution_clock::now() - progress_bar.start_time_;
		double estimated_remaining_runtime{(runtime_so_far.count() * ((100.0 / progress_in_percent) - 1.0)) / (60.0 * 60.0)};
		if (estimated_remaining_runtime >= 10) {
			num_digits++;
		}
		if (estimated_remaining_runtime >= 100) {
			num_digits++;
		}
		if (estimated_remaining_runtime >= 1000) {
			num_digits++;
		}
		while(num_digits++ < max_num_digits) {
			os << " ";
		}
		os << " " << estimated_remaining_runtime;
	}
	else {
		os << "       ?";
	}
	os << " h rem.";
	os.unsetf(std::ios::floatfield);
	os.precision(precision);
	return os;
}


}

#endif /* SRC_ENV_PROGRESS_BAR_IPP_ */
