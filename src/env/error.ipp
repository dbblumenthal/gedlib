/*!
 * @file error.ipp
 * @brief ged::Error class definition.
 */

#ifndef SRC_ENV_ERROR_IPP_
#define SRC_ENV_ERROR_IPP_

namespace ged {

Error ::
Error(const std::string & message) :
std::runtime_error(message) {}

}

#endif /* SRC_ENV_ERROR_IPP_ */
