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

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/tests/sspr2018/bin)

add_custom_target(sspr2018 DEPENDS learn_ring_params learn_walks_depths learn_subgraph_depths test_lsape_based_methods)

add_executable(learn_ring_params src/learn_ring_params.cpp)
target_link_libraries(learn_ring_params gxlgedlib)
if(APPLE)
  add_custom_command(TARGET learn_ring_params POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/learn_ring_params)
  add_custom_command(TARGET learn_ring_params POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/learn_ring_params)
  add_custom_command(TARGET learn_ring_params POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/learn_ring_params)
endif()

add_executable(learn_walks_depths src/learn_walks_depths.cpp)
target_link_libraries(learn_walks_depths gxlgedlib)
if(APPLE)
  add_custom_command(TARGET learn_walks_depths POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/learn_walks_depths)
  add_custom_command(TARGET learn_walks_depths POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/learn_walks_depths)
  add_custom_command(TARGET learn_walks_depths POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/learn_walks_depths)
endif()

add_executable(learn_subgraph_depths src/learn_subgraph_depths.cpp)
target_link_libraries(learn_subgraph_depths gxlgedlib)
if(APPLE)
  add_custom_command(TARGET learn_subgraph_depths POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/learn_subgraph_depths)
  add_custom_command(TARGET learn_subgraph_depths POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/learn_subgraph_depths)
  add_custom_command(TARGET learn_subgraph_depths POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/learn_subgraph_depths)
endif()

add_executable(test_lsape_based_methods src/test_lsape_based_methods.cpp)
target_link_libraries(test_lsape_based_methods gxlgedlib)
if(APPLE)
  add_custom_command(TARGET test_lsape_based_methods POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/test_lsape_based_methods)
  add_custom_command(TARGET test_lsape_based_methods POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/test_lsape_based_methods)
  add_custom_command(TARGET test_lsape_based_methods POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/test_lsape_based_methods)
endif()

