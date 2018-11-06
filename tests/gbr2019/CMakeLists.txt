############################################################################
#                                                                          #
#   Copyright (C) 2018 by David B. Blumenthal                              #
#                                                                          #
#   This file is part of GEDLIB.                                           #
#                                                                          #
#   GEDLIB is free software: you can redistribute it and/or modify it      #
#   under the terms of the GNU Lesser General Public License as published  #
#   by the Free Software Foundation, either version 3 of the License, or   #
#   (at your option) any later version.                                    #
#                                                                          #
#   GEDLIB is distributed in the hope that it will be useful,              #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           #
#   GNU Lesser General Public License for more details.                    #
#                                                                          #
#   You should have received a copy of the GNU Lesser General Public       #
#   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. #
#                                                                          #
############################################################################

# Set executable path.

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/tests/gbr2019/bin)

# Add targets.

add_custom_target(gbr2019 DEPENDS test_fgw_ged)

add_executable(test_fgw_ged src/test_fgw_ged.cpp)
target_link_libraries(test_fgw_ged gxlgedlib)
target_link_libraries(test_fgw_ged fgwgedlib)
if(APPLE)
  add_custom_command(TARGET test_fgw_ged POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/test_fgw_ged)
  add_custom_command(TARGET test_fgw_ged POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/test_fgw_ged)
  add_custom_command(TARGET test_fgw_ged POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/test_fgw_ged)
  add_custom_command(TARGET test_fgw_ged POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/test_fgw_ged)
  add_custom_command(TARGET test_fgw_ged POST_BUILD COMMAND install_name_tool -change libgurobi80.so ${GUROBI_HOME}/lib/libgurobi80.so ${EXECUTABLE_OUTPUT_PATH}/test_fgw_ged)
endif()