cmake_minimum_required(VERSION 3.27)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(mpm C CXX)
set(CMAKE_CXX_STANDARD 17)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Options
option(MPM_TESTING_ONLY "Only compile lib and tests" OFF)

# Libigl also includes Eigen, etc
include(libigl)

# OpenMP
# (TO-DO: replace with TBB)
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions(-DOMP_NESTED)
endif()

# Simulation library
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

# Interactive viewer
if (NOT MPM_TESTING_ONLY)
	igl_include(glfw)
	add_executable(sphere ${CMAKE_CURRENT_SOURCE_DIR}/test/sphere.cpp)
	target_link_libraries(sphere PUBLIC mpm_optimization igl::glfw)
	target_include_directories(sphere PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()

# Unit testing (TO-DO)
enable_testing()
add_executable(test_solver ${CMAKE_CURRENT_SOURCE_DIR}/test/solver.cpp)
target_link_libraries(test_solver PUBLIC mpm_optimization)
target_include_directories(test_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_test(NAME TestSolver COMMAND test_solver)