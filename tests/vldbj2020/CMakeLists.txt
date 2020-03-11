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

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/tests/vldbj2020/bin)

# Add targets.

add_custom_target(vldbj2020 DEPENDS vldbj_train_ring vldbj_train_walks vldbj_train_subgraph vldbj_train_ml vldbj_test_lsape_based_methods vldbj_test_best_methods)

add_executable(vldbj_train_ring src/train_ring.cpp)
target_link_libraries(vldbj_train_ring gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_train_ring POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ring)
  add_custom_command(TARGET vldbj_train_ring POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ring)
  add_custom_command(TARGET vldbj_train_ring POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ring)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_train_ring POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ring)
    add_custom_command(TARGET vldbj_train_ring POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ring)
  endif()
endif()

add_executable(vldbj_train_walks src/train_walks.cpp)
target_link_libraries(vldbj_train_walks gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_train_walks POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_walks)
  add_custom_command(TARGET vldbj_train_walks POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_walks)
  add_custom_command(TARGET vldbj_train_walks POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_walks)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_train_walks POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_walks)
    add_custom_command(TARGET vldbj_train_walks POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_walks)
  endif()
endif()

add_executable(vldbj_train_subgraph src/train_subgraph.cpp)
target_link_libraries(vldbj_train_subgraph gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_train_subgraph POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_subgraph)
  add_custom_command(TARGET vldbj_train_subgraph POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_subgraph)
  add_custom_command(TARGET vldbj_train_subgraph POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_subgraph)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_train_subgraph POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_subgraph)
    add_custom_command(TARGET vldbj_train_subgraph POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_subgraph)
  endif()
endif()

add_executable(vldbj_train_ml src/train_ml.cpp)
target_link_libraries(vldbj_train_ml gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_train_ml POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ml)
  add_custom_command(TARGET vldbj_train_ml POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ml)
  add_custom_command(TARGET vldbj_train_ml POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ml)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_train_ml POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ml)
    add_custom_command(TARGET vldbj_train_ml POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_train_ml)
  endif()
endif()

add_executable(vldbj_test_lsape_based_methods src/test_lsape_based_methods.cpp)
target_link_libraries(vldbj_test_lsape_based_methods gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_test_lsape_based_methods POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lsape_based_methods)
  add_custom_command(TARGET vldbj_test_lsape_based_methods POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lsape_based_methods)
  add_custom_command(TARGET vldbj_test_lsape_based_methods POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lsape_based_methods)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_test_lsape_based_methods POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lsape_based_methods)
    add_custom_command(TARGET vldbj_test_lsape_based_methods POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lsape_based_methods)
  endif()
endif()

add_executable(vldbj_test_lp_based_methods src/test_lp_based_methods.cpp)
target_link_libraries(vldbj_test_lp_based_methods gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_test_lp_based_methods POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lp_based_methods)
  add_custom_command(TARGET vldbj_test_lp_based_methods POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lp_based_methods)
  add_custom_command(TARGET vldbj_test_lp_based_methods POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lp_based_methods)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_test_lp_based_methods POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lp_based_methods)
    add_custom_command(TARGET vldbj_test_lp_based_methods POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_lp_based_methods)
  endif()
endif()

add_executable(vldbj_test_best_methods src/test_best_methods.cpp)
target_link_libraries(vldbj_test_best_methods gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_test_best_methods POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_best_methods)
  add_custom_command(TARGET vldbj_test_best_methods POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_best_methods)
  add_custom_command(TARGET vldbj_test_best_methods POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_best_methods)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_test_best_methods POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_best_methods)
    add_custom_command(TARGET vldbj_test_best_methods POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_best_methods)
  endif()
endif()

add_executable(vldbj_test_ls_based_methods src/test_ls_based_methods.cpp)
target_link_libraries(vldbj_test_ls_based_methods gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_test_ls_based_methods POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_ls_based_methods)
  add_custom_command(TARGET vldbj_test_ls_based_methods POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_ls_based_methods)
  add_custom_command(TARGET vldbj_test_ls_based_methods POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_ls_based_methods)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_test_ls_based_methods POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_ls_based_methods)
    add_custom_command(TARGET vldbj_test_ls_based_methods POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_ls_based_methods)
  endif()
endif()

add_executable(vldbj_test_misc_methods src/test_misc_methods.cpp)
target_link_libraries(vldbj_test_misc_methods gxlgedlib)
if(APPLE)
  add_custom_command(TARGET vldbj_test_misc_methods POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_misc_methods)
  add_custom_command(TARGET vldbj_test_misc_methods POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_misc_methods)
  add_custom_command(TARGET vldbj_test_misc_methods POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_misc_methods)
  if(GUROBI_HOME)
    add_custom_command(TARGET vldbj_test_misc_methods POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_misc_methods)
    add_custom_command(TARGET vldbj_test_misc_methods POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/vldbj_test_misc_methods)
  endif()
endif()