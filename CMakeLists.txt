# Copyright (c) 2016 University of Minnesota
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
#
# By Matt Overby (http://www.mattoverby.net)

cmake_minimum_required(VERSION 2.6)
project(mpm C CXX)

#
#	Build Mode
#
set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
set( CMAKE_BUILD_TYPE Release )

#
#	Variables
#

add_definitions( -DMPM_SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}" )


############################################################
#
#	Libraries
#
############################################################


# Eigen3
find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions( -DOMP_NESTED )
endif()

# CppOptimizationLibrary (header only)
set( CPPMIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/cppoptlib" )
include_directories( ${CPPMIN_DIR}/include )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )
set( MPM_SRCS
	src/Solver.hpp			src/Solver.cpp
	src/MPM.hpp			src/MPM.cpp
	src/Interp.hpp
	src/Particle.hpp
)


############################################################
#
#	Binaries
#
############################################################

add_executable( sphere src/sphere.cpp ${MPM_SRCS} )
target_link_libraries( sphere ${MPM_LIBRARIES} )


