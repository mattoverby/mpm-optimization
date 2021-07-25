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

############################################################
#
#	Libraries
#
############################################################

# libigl for rendering
include(DownloadProject)
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF) # off = header only
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         OFF)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
download_project(PROJ libigl
	GIT_REPOSITORY https://github.com/libigl/libigl.git
	GIT_TAG main
	UPDATE_DISCONNECTED 1
	QUIET)
add_subdirectory(${libigl_SOURCE_DIR} ${libigl_BINARY_DIR})

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions(-DOMP_NESTED)
endif()

# Include headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/deps/mcloptlib/include)
include_directories(SYSTEM ${glad_SOURCE_DIR}/include)
include_directories(SYSTEM ${libigl_SOURCE_DIR}/include)
include_directories(SYSTEM ${libigl_SOURCE_DIR}/external/eigen)

set(MPM_SRCS
	src/Solver.hpp
	src/Solver.cpp
	src/MPM.hpp
	src/MPM.cpp
	src/Interp.hpp
	src/Particle.hpp)

############################################################
#
#	Binaries
#
############################################################

add_executable(sphere src/sphere.cpp ${MPM_SRCS})
target_link_libraries(sphere glad glfw)
