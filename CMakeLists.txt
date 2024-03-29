#
# Copyright (C) 2015   Michael Lahnert
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

project (NUMSIM)

cmake_minimum_required (VERSION 2.8)

# force to set a build type
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# set compiler
set (CMAKE_CC_COMPILER "/usr/bin/gcc")
set (CMAKE_CXX_COMPILER "/usr/bin/g++")

# set some good compiler flags
#set (GOOD_FLAGS "-Wall -Wextra -pedantic -g -O3 -fno-omit-frame-pointer -march=native -std=c++11")
set (GOOD_FLAGS "-g -std=c++11")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOOD_FLAGS}")

# debug visu?
option (DEBUG_VISU
  "Enable debugging visualization" ON)

# find sdl 2
#INCLUDE(FindPkgConfig)
#PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
#INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${SDL2_INCLUDE_DIRS})

# find mpi
find_package (MPI REQUIRED)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${MPI_INCLUDE_PATH})

find_package(CUDA QUIET REQUIRED)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS})
list(APPEND CUDA_NVCC_FLAGS "-g -O3")

# Using glob is discouraged in cmake.
set (HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/src/communicator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/typedef.hpp
  #${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/particle.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/falggrid.hpp
  )

set (SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/communicator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.cpp
  #${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/particle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/flaggrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mg_impl.cu
  )

# add header and source of visu if it shall be built
# additionally set a flag
#if (DEBUG_VISU)
#  set (HEADERS
#    ${HEADERS}
#    ${CMAKE_CURRENT_SOURCE_DIR}/src/visu.hpp
#    )
#
#  set (SOURCES
#    ${SOURCES}
#    ${CMAKE_CURRENT_SOURCE_DIR}/src/visu.cpp
#    )
#
#  add_definitions(-DUSE_DEBUG_VISU)
#endif (DEBUG_VISU)

# build executable and link libs against it
cuda_add_executable(numsim ${SOURCES})
target_link_libraries (numsim ${SDL2_LIBRARIES} ${MPI_LIBRARIES} cuda)

if (MPI_COMPILE_FLAGS)
  set_target_properties(numsim PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if (MPI_LINK_FLAGS)
  set_target_properties (numsim PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
