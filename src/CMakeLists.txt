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

add_library(gxlgedlib SHARED env/ged_env.gxl.cpp)
set_target_properties(gxlgedlib PROPERTIES SUFFIX ".so")
target_link_libraries(gxlgedlib nomad doublefann svm)
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
