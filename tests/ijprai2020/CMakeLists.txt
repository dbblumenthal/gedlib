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

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/tests/ijprai2020/bin)

# Add targets.

add_custom_target(ijprai2020 DEPENDS train_ring train_walks train_subgraph train_ml run_tests run_ipfp)

add_executable(train_ring src/train_ring.cpp)
target_link_libraries(train_ring gxlgedlib)
if(APPLE)
  add_custom_command(TARGET train_ring POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/train_ring)
  add_custom_command(TARGET train_ring POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/train_ring)
  add_custom_command(TARGET train_ring POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/train_ring)
endif()

add_executable(train_walks src/train_walks.cpp)
target_link_libraries(train_walks gxlgedlib)
if(APPLE)
  add_custom_command(TARGET train_walks POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/train_walks)
  add_custom_command(TARGET train_walks POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/train_walks)
  add_custom_command(TARGET train_walks POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/train_walks)
endif()

add_executable(train_subgraph src/train_subgraph.cpp)
target_link_libraries(train_subgraph gxlgedlib)
if(APPLE)
  add_custom_command(TARGET train_subgraph POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/train_subgraph)
  add_custom_command(TARGET train_subgraph POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/train_subgraph)
  add_custom_command(TARGET train_subgraph POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/train_subgraph)
endif()

add_executable(train_ml src/train_ml.cpp)
target_link_libraries(train_ml gxlgedlib)
if(APPLE)
  add_custom_command(TARGET train_ml POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/train_ml)
  add_custom_command(TARGET train_ml POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/train_ml)
  add_custom_command(TARGET train_ml POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/train_ml)
endif()

add_executable(run_tests src/tests.cpp)
target_link_libraries(run_tests gxlgedlib)
if(APPLE)
  add_custom_command(TARGET run_tests POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/run_tests)
  add_custom_command(TARGET run_tests POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/run_tests)
  add_custom_command(TARGET run_tests POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/run_tests)
endif()

add_executable(run_ipfp src/run_ipfp.cpp)
target_link_libraries(run_ipfp gxlgedlib)
if(APPLE)
  add_custom_command(TARGET run_ipfp POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/run_ipfp)
  add_custom_command(TARGET run_ipfp POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/run_ipfp)
  add_custom_command(TARGET run_ipfp POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/run_ipfp)
endif()
