# FindOpenMPCustom.cmake
# Custom OpenMP finder for macOS with Homebrew libomp support
#
# This module defines:
#   OpenMP_FOUND - True if OpenMP is found
#   OpenMP_CXX_FOUND - True if OpenMP for CXX is found
#   OpenMP::OpenMP_CXX - Imported target for OpenMP CXX support
#
# Usage:
#   find_package(OpenMPCustom)
#   if(OpenMP_CXX_FOUND)
#       target_link_libraries(your_target PRIVATE OpenMP::OpenMP_CXX)
#   endif()

# First try the standard FindOpenMP
find_package(OpenMP QUIET)

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found using standard CMake module")
    set(OpenMPCustom_FOUND TRUE)
else()
    message(STATUS "Standard OpenMP not found, trying to find libomp (Homebrew)")
    
    # Try to find libomp (Homebrew version)
    find_path(OpenMP_INCLUDE_DIR omp.h 
        HINTS 
        /opt/homebrew/Cellar/libomp/21.1.1/include
        /opt/homebrew/include
        /usr/local/include
        PATHS
        /opt/homebrew/Cellar/libomp/*/include
        /usr/local/Cellar/libomp/*/include
        PATH_SUFFIXES include
    )
    
    find_library(OpenMP_LIBRARY 
        NAMES omp gomp iomp5
        HINTS 
        /opt/homebrew/Cellar/libomp/21.1.1/lib
        /opt/homebrew/lib
        /usr/local/lib
        PATHS
        /opt/homebrew/Cellar/libomp/*/lib
        /usr/local/Cellar/libomp/*/lib
        PATH_SUFFIXES lib
    )
    
    if(OpenMP_INCLUDE_DIR AND OpenMP_LIBRARY)
        set(OpenMP_FOUND TRUE)
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIBRARIES ${OpenMP_LIBRARY})
        message(STATUS "Found OpenMP (libomp): ${OpenMP_LIBRARY}")
        
        # Create imported target if it doesn't already exist
        if(NOT TARGET OpenMP::OpenMP_CXX)
            add_library(OpenMP::OpenMP_CXX SHARED IMPORTED)
            set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                IMPORTED_LOCATION ${OpenMP_LIBRARY}
                INTERFACE_INCLUDE_DIRECTORIES ${OpenMP_INCLUDE_DIR}
                INTERFACE_COMPILE_OPTIONS "-Xpreprocessor;-fopenmp"
            )
        endif()
        set(OpenMP_CXX_FOUND TRUE)
        set(OpenMPCustom_FOUND TRUE)
    else()
        message(WARNING "OpenMP not found. Parallel code will not be optimized.")
        set(OpenMP_FOUND FALSE)
        set(OpenMP_CXX_FOUND FALSE)
        set(OpenMPCustom_FOUND FALSE)
    endif()
endif()

# Handle standard CMake find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMPCustom
    FOUND_VAR OpenMPCustom_FOUND
    REQUIRED_VARS OpenMP_CXX_FOUND
)
