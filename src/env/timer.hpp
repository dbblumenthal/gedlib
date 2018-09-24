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
