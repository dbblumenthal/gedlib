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
