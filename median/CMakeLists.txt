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

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/median/bin)

# Add targets.

add_custom_target(median DEPENDS median_letter median_chem median_tests)

add_custom_target(cluster DEPENDS cluster_letter)

add_custom_target(bst DEPENDS bst_init_letter bst_query_letter)

add_executable(median_letter tests/median_letter.cpp)
target_link_libraries(median_letter gxlgedlib)
if(APPLE)
  add_custom_command(TARGET median_letter POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/median_letter)
  add_custom_command(TARGET median_letter POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/median_letter)
  add_custom_command(TARGET median_letter POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/median_letter)
endif()

add_executable(median_chem tests/median_chem.cpp)
target_link_libraries(median_chem gxlgedlib)
if(APPLE)
  add_custom_command(TARGET median_chem POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/median_chem)
  add_custom_command(TARGET median_chem POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/median_chem)
  add_custom_command(TARGET median_chem POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/median_chem)
endif()

add_executable(median_tests tests/median_tests.cpp)
target_link_libraries(median_tests gxlgedlib)
if(APPLE)
  add_custom_command(TARGET median_tests POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/median_tests)
  add_custom_command(TARGET median_tests POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/median_tests)
  add_custom_command(TARGET median_tests POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/median_tests)
endif()

add_executable(median_letter_demo tests/median_letter_demo.cpp)
target_link_libraries(median_letter_demo gxlgedlib)
if(APPLE)
  add_custom_command(TARGET median_letter_demo POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/median_letter_demo)
  add_custom_command(TARGET median_letter_demo POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/median_letter_demo)
  add_custom_command(TARGET median_letter_demo POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/median_letter_demo)
endif()

add_executable(cluster_letter tests/cluster_letter.cpp)
target_link_libraries(cluster_letter gxlgedlib)
if(APPLE)
  add_custom_command(TARGET cluster_letter POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/cluster_letter)
  add_custom_command(TARGET cluster_letter POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/cluster_letter)
  add_custom_command(TARGET cluster_letter POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/cluster_letter)
endif()

add_executable(bst_init_letter tests/bst_init_letter.cpp)
target_link_libraries(bst_init_letter gxlgedlib)
if(APPLE)
  add_custom_command(TARGET bst_init_letter POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/bst_init_letter)
  add_custom_command(TARGET bst_init_letter POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/bst_init_letter)
  add_custom_command(TARGET bst_init_letter POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/bst_init_letter)
endif()

add_executable(bst_query_letter tests/bst_query_letter.cpp)
target_link_libraries(bst_query_letter gxlgedlib)
if(APPLE)
  add_custom_command(TARGET bst_query_letter POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/bst_query_letter)
  add_custom_command(TARGET bst_query_letter POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/bst_query_letter)
  add_custom_command(TARGET bst_query_letter POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/bst_query_letter)
endif()

