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
 * @file ged_env.gxl.cpp
 * @brief ged::GEDEnv<ged::GXLUserNodeID, ged::GXLLabel, ged::GXLLabel> template instantiation.
 */

#ifndef SRC_ENV_GED_ENV_GXL_CPP_
#define SRC_ENV_GED_ENV_GXL_CPP_

#include "ged_env.hpp"

namespace ged {

template class GEDEnv<GXLNodeID, GXLLabel, GXLLabel>;

}

#endif /* SRC_ENV_GED_ENV_GXL_CPP_ */
