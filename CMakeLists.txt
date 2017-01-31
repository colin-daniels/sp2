cmake_minimum_required(VERSION 3.2)
project(sp2 CXX)

# support unix only
if(NOT UNIX)
    message(FATAL_ERROR "SP2 is supported on unix only." )
endif()

# set default build type to release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: \n\
None Debug Release RelWithDebInfo MinSizeRel Coverage."
            FORCE)
endif()

# add standard modules
include(ExternalProject)
include(CMakeParseArguments)
include(ProcessorCount)

# as well as our custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
include(PreventInSourceBuilds)
include(DownloadProject)
include(FindPythonModule)

# enable ctest
enable_testing()

################################################################################
## Options                                                                    ##
################################################################################

option(SP2_ENABLE_MPI       "Enable MPI" ON)
option(SP2_ENABLE_LAMMPS    "Enable LAMMPS" ON)
option(SP2_ENABLE_PHONOPY   "Enable Phonopy" ON)
option(SP2_ENABLE_TESTS     "Enable/build tests" ON)

option(BUILD_SHARED_LIBS    "Build shared libraries (DLLs)." OFF)

################################################################################
## Source files                                                               ##
################################################################################

set(SOURCE_FILES
        src/common/vec3_t.hpp
        src/common/vec3_t.cpp
        src/common/function_types.hpp
        src/common/function_types.cpp
        src/common/structure_t.cpp
        src/common/structure_t.hpp
        src/common/util/rotations.cpp
        src/common/util/rotations.hpp
        src/common/util/blas.hpp
        src/common/util/random.cpp
        src/common/util/random.hpp
        src/common/util/templates.hpp
        src/common/util/numerical_diff.hpp
        src/common/util/interpolate.cpp
        src/common/util/math.hpp
        src/common/util/timing.hpp
        src/common/util/timing.cpp
        src/common/io/structure.cpp
        src/common/io/structure.hpp
        src/common/io/util.cpp
        src/common/io/util.hpp
        src/common/io/file_types/xyz.cpp
        src/common/io/file_types/xyz.hpp
        src/common/io/file_types/cif.cpp
        src/common/io/file_types/cif.hpp
        src/common/io/file_types/poscar.cpp
        src/common/io/file_types/poscar.hpp
        src/common/json/json.cpp
        src/common/json/json.hpp
        src/common/json/json_serializable_t.hpp
        src/common/minimize/minimize.hpp
        src/common/minimize/settings.cpp
        src/common/minimize/settings.hpp
        src/common/minimize/test/test_functions.cpp
        src/common/minimize/test/test_functions.hpp
        src/common/minimize/linesearch.cpp
        src/common/minimize/linear_cg.cpp
        src/common/minimize/acgsd.cpp
################################################################################
        IF_MPI_ENABLED src/common/util/mpi.hpp
        # PSO requires MPI
        IF_MPI_ENABLED src/common/minimize/pso/adaptive_pso.cpp
        IF_MPI_ENABLED src/common/minimize/pso/particle_t.cpp
        IF_MPI_ENABLED src/common/minimize/pso/particle_t.hpp
        IF_MPI_ENABLED src/common/minimize/pso/swarm_t.cpp
        IF_MPI_ENABLED src/common/minimize/pso/swarm_t.hpp
        # and by extension so does symmetry search
        IF_MPI_ENABLED src/run/run_symm.cpp
        IF_MPI_ENABLED src/symm/pso_adapters.cpp
        IF_MPI_ENABLED src/symm/pso_adapters.hpp
        IF_MPI_ENABLED src/symm/util.cpp
        IF_MPI_ENABLED src/symm/util.hpp
        IF_MPI_ENABLED src/symm/space_group_t.cpp
        IF_MPI_ENABLED src/symm/space_group_t.hpp
        IF_MPI_ENABLED src/symm/symm_settings_t.cpp
        IF_MPI_ENABLED src/symm/symm_settings_t.hpp
        IF_MPI_ENABLED src/symm/system_control_t.cpp
        IF_MPI_ENABLED src/symm/system_control_t.hpp
################################################################################
        IF_LAMMPS_ENABLED src/lammps/system_control_t.cpp
        IF_LAMMPS_ENABLED src/lammps/lammps_interface.hpp
        IF_LAMMPS_ENABLED src/lammps/settings_t.cpp
        IF_LAMMPS_ENABLED src/lammps/settings_t.hpp
################################################################################
        IF_PHONOPY_ENABLED src/run/run_phonopy.cpp
        IF_PHONOPY_ENABLED src/phonopy/phonopy_settings.cpp
        IF_PHONOPY_ENABLED src/phonopy/phonopy_settings.hpp
################################################################################
        src/common/graph/ud_graph_t.cpp
        src/common/graph/ud_graph_t.hpp
        src/common/graph/permutation_t.cpp
        src/common/graph/permutation_t.hpp
        src/common/graph/graph.cpp
        src/common/graph/graph.hpp
        src/common/graph/test/test_ud_graph_t.cpp
        src/common/neighbor/utility_functions.cpp
        src/common/neighbor/utility_functions.hpp
        src/common/neighbor/periodic_cell_t.hpp
        src/common/neighbor/cell_array_t.cpp
        src/common/neighbor/cell_array_t.hpp
        src/common/neighbor/bond_control_t.cpp
        src/common/neighbor/bond_control_t.hpp
        src/common/enums.cpp
        src/common/enums.hpp
        src/common/neighbor/neighbor_control_t.hpp
        src/common/neighbor/neighbor_control_t.cpp
        src/airebo/system_control_t.cpp
        src/airebo/system_control_t.hpp
        src/airebo/interpolation_coeff.hpp
        src/airebo/utility_functions.hpp
        src/airebo/test.hpp
        src/airebo/test.cpp
        # src/atac/system_control_t.cpp
        # src/atac/system_control_t.hpp
        src/atac/settings_t.hpp
        src/atac/settings_t.cpp
        src/run/common.cpp
        # src/run/run_relaxation.cpp
        # src/run/run_atac.cpp
        src/run/run_types.hpp
        src/run/run_settings_t.cpp
        src/run/run_settings_t.hpp
        src/common/json/json_serializable_t.cpp
        src/airebo/airebo_util.hpp)

include(Sp2ParseSources)
parse_sources(SOURCE_FILES ${SOURCE_FILES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR})

# add subdirectories for external code as well as the code generation target
add_subdirectory(external)
add_subdirectory(src/common/codegen)

link_libraries(templ)
link_libraries(jsoncpp)

if(SP2_ENABLE_LAMMPS)
    link_libraries(lammps)
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_definitions(-DSP2_DEBUG)
endif()

if(SP2_ENABLE_TESTS)
    add_test(all sp2 --test)
    add_definitions(-DSP2_ENABLE_TESTS)
    link_libraries(gtest)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

################################################################################
## MPI                                                                        ##
################################################################################

if(SP2_ENABLE_MPI)
    find_package(MPI REQUIRED)

    add_definitions(-DSP2_ENABLE_MPI)
    list(APPEND CMAKE_CXX_COMPILE_FLAGS
            ${MPI_CXX_COMPILE_FLAGS})

    include_directories(${MPI_CXX_INCLUDE_PATH})

    set(CMAKE_CXX_LINK_FLAGS
            ${CMAKE_CXX_LINK_FLAGS}
            ${MPI_CXX_LINK_FLAGS})

    link_libraries(boost_mpi
            boost_serialization
            ${MPI_CXX_LIBRARIES})
else()
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)

    link_libraries(Threads::Threads)
endif()

################################################################################
## Phonopy                                                                    ##
################################################################################

function(find_phonopy version)
    cmake_parse_arguments(FP_ARGS "REQUIRED" "" "" ${ARGN})

    # force find_package() to look for python again
    unset(PYTHON_EXECUTABLE CACHE)
    unset(PYTHONINTERP_FOUND CACHE)

    if(FP_ARGS_REQUIRED)
        set(REQUIRED_STR "REQUIRED")
    endif()

    find_package(PythonInterp ${version} ${REQUIRED_STR})
    if(PYTHONINTERP_FOUND)
        find_python_module(phonopy ${REQUIRED_STR})
    endif()
endfunction(find_phonopy)

if(SP2_ENABLE_PHONOPY AND NOT PYMOD_PHONOPY_LOCATION)
    # first we try to find the python 3+ version of phonopy
    find_phonopy(3)
    # then if that fails, fall back to 2.7
    if(NOT PYMOD_PHONOPY_FOUND)
        find_phonopy(2.7)
    endif()

    add_definitions(-DSP2_ENABLE_PHONOPY)
endif()

################################################################################
## Main targets                                                               ##
################################################################################

# sp2 executable (sp2)
add_executable(sp2-bin
        src/main.cpp
        ${SOURCE_FILES})

# explicitly set executable filename to 'sp2'
set_target_properties(sp2-bin PROPERTIES
        OUTPUT_NAME sp2)

# sp2 library (libsp2.a/so)
add_library(sp2-lib
        src/run/lib/library_interface.cpp
        ${SOURCE_FILES})

# explicitly set library filename to libsp2.so/libsp2.a
# additionally, stop the library from being made by default
set_target_properties(sp2-lib PROPERTIES
        OUTPUT_NAME sp2
        EXCLUDE_FROM_ALL 1
        EXCLUDE_FROM_DEFAULT_BUILD 1)

# code generation dependencies
add_dependencies(sp2-bin code_generation)
add_dependencies(sp2-lib code_generation)

