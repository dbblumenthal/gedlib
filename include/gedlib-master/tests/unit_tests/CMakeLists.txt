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

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/tests/unit_tests/bin)

add_custom_target(unit_tests DEPENDS ged_env_tests lsap_solver_tests)

add_executable(ged_env_tests src/catch.cpp src/ged_env_test.cpp)
target_link_libraries(ged_env_tests gxlgedlib)
if(APPLE)
  add_custom_command(TARGET ged_env_tests POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/ged_env_tests)
  add_custom_command(TARGET ged_env_tests POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/ged_env_tests)
  add_custom_command(TARGET ged_env_tests POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/ged_env_tests)
endif()

add_executable(lsap_solver_tests src/catch.cpp src/lsap_solver_test.cpp)