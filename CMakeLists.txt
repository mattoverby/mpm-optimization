cmake_minimum_required(VERSION 3.27)

# Generate compile_commands.json for language servers and IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(mpm LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# Add cmake modules to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ###########################################################
#
# Libraries
#
# ###########################################################

# Find and link dependencies using Vcpkg
find_package(glad CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
# libigl[core,opengl,glfw]
find_package(libigl CONFIG REQUIRED)

# Find OpenMP using custom module
find_package(OpenMPCustom)

# ###########################################################
#
# Binaries
#
# ###########################################################

add_executable(sphere)

file(GLOB_RECURSE MPM_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp")
target_sources(sphere PRIVATE ${MPM_SRCS})

# Link libraries
target_link_libraries(sphere PRIVATE glad::glad glfw)
target_link_libraries(sphere PRIVATE igl::igl_core igl::igl_glfw igl::igl_opengl igl_copyleft::igl_copyleft_core)
target_link_libraries(sphere PRIVATE Eigen3::Eigen)

# Link OpenMP if found
if(OpenMP_CXX_FOUND)
    target_link_libraries(sphere PRIVATE OpenMP::OpenMP_CXX)
endif()

# Include directories
target_include_directories(sphere
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/deps/mcloptlib/include
)

target_compile_definitions(sphere PRIVATE MPM_SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

# IDE support
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src PREFIX "Source Files" FILES ${MPM_SRCS})