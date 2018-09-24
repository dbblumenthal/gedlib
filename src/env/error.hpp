/*!
 * @file error.hpp
 * @brief ged::Error class declaration.
 */

#ifndef SRC_ENV_ERROR_HPP_
#define SRC_ENV_ERROR_HPP_

#include <stdexcept>

namespace ged {

/*!
 * @brief Runtime error class.
 */
class Error : public std::runtime_error {

public:

	/*!
	 * @brief Constructor.
	 * @param[in] message Error message.
	 */
	Error(const std::string & message);

};

}

#include "error.ipp"

#endif /* SRC_ENV_ERROR_HPP_ */
