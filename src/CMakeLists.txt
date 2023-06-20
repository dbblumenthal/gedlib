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

set(CPM_DOWNLOAD_VERSION 0.34.0)

if(CPM_SOURCE_CACHE)
  set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
elseif(DEFINED ENV{CPM_SOURCE_CACHE})
  set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
else()
  set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
endif()

add_library(gxlgedlib SHARED env/ged_env.gxl.cpp)
set_target_properties(gxlgedlib PROPERTIES SUFFIX ".so")
target_link_libraries(gxlgedlib nomad doublefann svm)

add_library(cuda_acc STATIC methods/branch_fast.cu)
target_compile_definitions(cuda_acc PUBLIC -D CUDA_FORCE_CDP1_IF_SUPPORTED)
add_compile_options(-fPIC -fno-pie)
set_target_properties(cuda_acc PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(gxlgedlib cuda_acc)

if(GUROBI_HOME)
  target_link_libraries(gxlgedlib ${GUROBI_STATLIB} ${GUROBI_DYLIB})
endif()

if(APPLE)
  add_custom_command(TARGET gxlgedlib POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${LIBRARY_OUTPUT_PATH}/libgxlgedlib.so)
  add_custom_command(TARGET gxlgedlib POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${LIBRARY_OUTPUT_PATH}/libgxlgedlib.so)
  add_custom_command(TARGET gxlgedlib POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${LIBRARY_OUTPUT_PATH}/libgxlgedlib.so)
  if(GUROBI_HOME)
    add_custom_command(TARGET gxlgedlib POST_BUILD COMMAND install_name_tool -change lib${GUROBI_STATLIB}.a ${GUROBI_HOME}/lib/lib${GUROBI_STATLIB}.a ${LIBRARY_OUTPUT_PATH}/libgxlgedlib.so)
    add_custom_command(TARGET gxlgedlib POST_BUILD COMMAND install_name_tool -change lib${GUROBI_DYLIB}.dylib ${GUROBI_HOME}/lib/lib${GUROBI_DYLIB}.dylib ${LIBRARY_OUTPUT_PATH}/libgxlgedlib.so)
  endif()
  target_link_libraries(gxlgedlib omp)
  add_custom_command(TARGET gxlgedlib POST_BUILD COMMAND install_name_tool -change libomp.dylib ${OMP_HOME}/lib/libomp.dylib ${LIBRARY_OUTPUT_PATH}/libgxlgedlib.so)
endif()

if(GED_ENABLE_CUDA)

endif()
