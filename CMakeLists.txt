cmake_minimum_required(VERSION 3.2)
project(sp2 CXX)
set(${PROJECT_NAME}_VERSION 0.1)

# add our custom module directory to the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

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

# interface target to carry along compile definitions related to project options
add_library(project-options INTERFACE)
# for each option, if it is enabled, it is added to the project-options target
# as an interface compile definition
include(VarsToDefines)
vars_to_defines(project-options INTERFACE
        SP2_ENABLE_MPI
        SP2_ENABLE_LAMMPS
        SP2_ENABLE_PHONOPY
        SP2_ENABLE_TESTS
        SP2_ENABLE_PYTHON)

# setup special per-file compilation flags
# This should be present on any .cpp file that imports CPython headers.
set(SP2_PYTHON_O_FLAGS -fno-strict-aliasing)

################################################################################
## Misc                                                                       ##
################################################################################

if(NOT UNIX)
    message(FATAL_ERROR "${PROJECT_NAME} is supported on unix only." )
endif()

# prevent people from building in-source (aka in the source tree)
include(PreventInSourceBuilds)

# set default build type to release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: \n\
None Debug Release RelWithDebInfo MinSizeRel Coverage."
            FORCE)
endif()

# enable ctest
enable_testing()

################################################################################
## Sources                                                                    ##
################################################################################

# documentation, target is 'docs', must be made explicitly (e.g. make docs)
add_subdirectory(docs EXCLUDE_FROM_ALL)

# interface target populated in the libs/ subdir in order to pass along new
# dependencies to the main target
add_library(project-deps INTERFACE)

# pass along project-options' compile definitions to the deps target
target_link_libraries(project-deps INTERFACE
        project-options)

# add subdirectories for external libraries/code, populates dependencies for
# the project-deps target
add_subdirectory(deps)

# main target directory
add_subdirectory(src)
