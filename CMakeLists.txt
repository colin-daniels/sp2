cmake_minimum_required(VERSION 3.2)
project(sp2 CXX)
set(${PROJECT_NAME}_VERSION 0.1)

################################################################################
## Options                                                                    ##
################################################################################

option(SP2_ENABLE_MPI                "Enable MPI" ON)
option(SP2_ENABLE_LAMMPS          "Enable LAMMPS" ON)
option(SP2_ENABLE_PHONOPY        "Enable Phonopy" ON)
option(SP2_ENABLE_TESTS      "Enable/build tests" OFF)
option(SP2_ENABLE_PYTHON "Enable Python bindings" ON)

option(SP2_BUILD_LIB     "Build library instead of binary" OFF)
option(BUILD_SHARED_LIBS  "Build shared libraries (DLLs)." OFF)

# setup special per-file compilation flags
# This should be present on any .cpp file that imports CPython headers.
set(SP2_PYTHON_O_FLAGS -fno-strict-aliasing)

if(NOT UNIX)
    message(FATAL_ERROR "${PROJECT_NAME} is supported on unix only." )
endif()

# set default build type to release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: \n\
None Debug Release RelWithDebInfo MinSizeRel Coverage."
            FORCE)
endif()

################################################################################
## Modules                                                                    ##
################################################################################

# add standard modules
include(CMakeParseArguments)
include(ProcessorCount)

# add custom module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(PreventInSourceBuilds)

# enable ctest
enable_testing()

################################################################################
## Sources                                                                    ##
################################################################################

# interface target populated in the libs/ subdir in order to pass along new
# dependencies to the main target
add_library(project-deps INTERFACE)
# add subdirectories for external libraries/code
add_subdirectory(deps)

# main target directory
add_subdirectory(src)
