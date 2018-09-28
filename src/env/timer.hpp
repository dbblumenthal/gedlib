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
 * @file timer.hpp
 * @brief ged::Timer class declaration.
 */

#ifndef SRC_ENV_TIMER_HPP_
#define SRC_ENV_TIMER_HPP_

#include <iostream>
#include <chrono>

namespace ged {

/*!
 * @brief A timer class that can be used by methods that support time limits.
 */
class Timer {
public:

	/*!
	 * @brief Constructs a timer for a given time limit.
	 * @param[in] time_limit_in_sec The time limit in seconds.
	 */
	Timer(double time_limit_in_sec);

	/*!
	 * @brief Checks if the time limit has expired.
	 * @return Boolean @p true if the time limit has expired and @p false otherwise.
	 */
	bool expired() const;

private:

	const double time_limit_in_sec_;

	const std::chrono::high_resolution_clock::time_point start_time_;
};

}

#include "timer.ipp"

#endif /* SRC_ENV_TIMER_HPP_ */
