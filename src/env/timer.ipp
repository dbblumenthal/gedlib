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
