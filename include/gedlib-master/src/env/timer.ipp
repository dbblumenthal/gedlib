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
 * @file timer.ipp
 * @brief ged::Timer class definition.
 */

#ifndef SRC_ENV_TIMER_IPP_
#define SRC_ENV_TIMER_IPP_

namespace ged {

Timer ::
Timer(double time_limit_in_sec) :
time_limit_in_sec_(time_limit_in_sec),
start_time_{std::chrono::high_resolution_clock::now()} {}

bool
Timer ::
expired() const {
	if (time_limit_in_sec_ > 0) {
		std::chrono::duration<double> runtime = std::chrono::high_resolution_clock::now() - start_time_;
		return (runtime.count() >= time_limit_in_sec_);
	}
	return false;
}

}

#endif /* SRC_ENV_TIMER_IPP_ */
