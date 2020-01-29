############################################################################
#                                                                          #
#   Copyright (C) 2020 by David B. Blumenthal                              #
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

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/app/IBD/bin)

add_custom_target(ibd DEPENDS compute_ibd_embeddings compute_ibd_prototypes tune_ibd_costs)

add_executable(compute_ibd_embeddings src/compute_ibd_embeddings.cpp)
target_link_libraries(compute_ibd_embeddings gxlgedlib)
if(APPLE)
  add_custom_command(TARGET compute_ibd_embeddings POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_embeddings)
  add_custom_command(TARGET compute_ibd_embeddings POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_embeddings)
  add_custom_command(TARGET compute_ibd_embeddings POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_embeddings)
  if(GUROBI_HOME)
    add_custom_command(TARGET compute_ibd_embeddings POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_embeddings)
    add_custom_command(TARGET compute_ibd_embeddings POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_embeddings)
  endif()
endif()

add_executable(compute_ibd_prototypes src/compute_ibd_prototypes.cpp)
target_link_libraries(compute_ibd_prototypes gxlgedlib)
if(APPLE)
  add_custom_command(TARGET compute_ibd_prototypes POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_prototypes)
  add_custom_command(TARGET compute_ibd_prototypes POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_prototypes)
  add_custom_command(TARGET compute_ibd_prototypes POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_prototypes)
  if(GUROBI_HOME)
    add_custom_command(TARGET compute_ibd_prototypes POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_prototypes)
    add_custom_command(TARGET compute_ibd_prototypes POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/compute_ibd_prototypes)
  endif()
endif()

add_executable(tune_ibd_costs src/tune_ibd_costs.cpp)
target_link_libraries(tune_ibd_costs gxlgedlib)
if(APPLE)
  add_custom_command(TARGET tune_ibd_costs POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/tune_ibd_costs)
  add_custom_command(TARGET tune_ibd_costs POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/tune_ibd_costs)
  add_custom_command(TARGET tune_ibd_costs POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/tune_ibd_costs)
  if(GUROBI_HOME)
    add_custom_command(TARGET tune_ibd_costs POST_BUILD COMMAND install_name_tool -change libgurobi_c++.a ${GUROBI_HOME}/lib/libgurobi_c++.a ${EXECUTABLE_OUTPUT_PATH}/tune_ibd_costs)
    add_custom_command(TARGET tune_ibd_costs POST_BUILD COMMAND install_name_tool -change lib${GUROBI_SHARED_LIB}.so ${GUROBI_HOME}/lib/lib${GUROBI_SHARED_LIB}.so ${EXECUTABLE_OUTPUT_PATH}/tune_ibd_costs)
  endif()
endif()