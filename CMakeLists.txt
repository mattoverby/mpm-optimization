# Copyright (c) 2016 Matt Overby
# 
# MPM-OPTIMIZATION Uses the BSD 2-Clause License (http://www.opensource.org/licenses/BSD-2-Clause)
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:  
# 1. Redistributions of source code must retain the above copyright notice, this list of
# conditions and the following disclaimer.  
# 2. Redistributions in binary form must reproduce the above copyright notice, this list
# of conditions and the following disclaimer in the documentation and/or other materials
# provided with the distribution.  
# THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF MINNESOTA, DULUTH OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.0)
project(mpm C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Release)

add_definitions(-DMPM_SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
option(MPM_TESTING_ONLY "Only compile lib and tests" OFF)

############################################################
#
#	Libraries
#
############################################################

# Libigl also includes Eigen, etc
include(libigl)

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions(-DOMP_NESTED)
endif()

set(MPM_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/src/LBFGS.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Solver.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Solver.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MPM.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MPM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Interp.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Particle.hpp)


add_library(mpm_optimization ${MPM_SRC})
target_link_libraries(mpm_optimization PUBLIC igl::core)

############################################################
#
#	Binaries
#
############################################################

if (NOT MPM_TESTING_ONLY)
	igl_include(glfw)
	add_executable(sphere ${CMAKE_CURRENT_SOURCE_DIR}/src/sphere.cpp)
	target_link_libraries(sphere PUBLIC mpm_optimization igl::glfw)
endif()

enable_testing()
add_executable(test_solver ${CMAKE_CURRENT_SOURCE_DIR}/test/solver.cpp)
target_link_libraries(test_solver PUBLIC mpm_optimization)
add_test(NAME TestSolver COMMAND test_solver)