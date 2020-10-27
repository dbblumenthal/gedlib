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

cmake_minimum_required(VERSION 2.6)

project(GEDLIB)
set(VERSION 1.0)

# Determine build type.

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Building GEDLIB with build type 'Release', as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release or Debug build type")
else()
  message(STATUS "Building GEDLIB with build type '${CMAKE_BUILD_TYPE}'.")
endif()

# Find Boost and Doxygen.

#find_package(Boost REQUIRED)
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doxyfile @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif()

# Set up the compiler.

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic -Wredundant-decls -Wdisabled-optimization -Wctor-dtor-privacy -Wnon-virtual-dtor -Werror=old-style-cast -Wsign-promo -Werror=return-type -Wno-unused-parameter")
if(APPLE)
  set(CMAKE_CXX_COMPILER clang++ CACHE STRING "C++ compiler" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()


set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_MACOSX_RPATH ON)

# Define variables.

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
set(LSAPE_HOME ${CMAKE_SOURCE_DIR}/ext/lsape.5)
set(NOMAD_HOME ${CMAKE_SOURCE_DIR}/ext/nomad.3.8.1)
set(LIBSVM_HOME ${CMAKE_SOURCE_DIR}/ext/libsvm.3.22)
set(EIGEN_HOME ${CMAKE_SOURCE_DIR}/ext/eigen.3.3.4)
set(FANN_HOME ${CMAKE_SOURCE_DIR}/ext/fann.2.2.0)
set(BOOST_HOME ${CMAKE_SOURCE_DIR}/ext/boost.1.69.0)

# Add include and link directories.

include_directories(SYSTEM ${BOOST_HOME} ${NOMAD_HOME}/ext/sgtelib/src ${NOMAD_HOME}/src ${LSAPE_HOME}/cpp/include ${EIGEN_HOME}/Eigen ${LIBSVM_HOME} ${FANN_HOME}/include ${GUROBI_HOME}/include ${OMP_HOME}/include)
link_directories(${NOMAD_HOME}/lib ${FANN_HOME}/lib ${LIBSVM_HOME} ${GUROBI_HOME}/lib ${OMP_HOME}/lib)

if(APPLE)
  include_directories(SYSTEM ${OMP_HOME}/include)
  link_directories(${OMP_HOME}/lib)
endif()

if(GUROBI_ROOT)
  add_definitions(-DGUROBI)
  if(APPLE)
    set(GUROBI_HOME ${GUROBI_ROOT}/mac64)
  else()
    set(GUROBI_HOME ${GUROBI_ROOT}/linux64)
  endif()
  include_directories(SYSTEM ${GUROBI_HOME}/include)
  link_directories(${GUROBI_HOME}/lib)
 endif()

# Add subdirectories.

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(median)