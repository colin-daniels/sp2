cmake_minimum_required(VERSION 3.2)
project(sp2 CXX)

################################################################################
## Options + misc                                                             ##
################################################################################

option(SP2_ENABLE_MPI                "Enable MPI" ON)
option(SP2_ENABLE_LAMMPS          "Enable LAMMPS" ON)
option(SP2_ENABLE_PHONOPY        "Enable Phonopy" ON)
option(SP2_ENABLE_TESTS      "Enable/build tests" OFF)
option(SP2_ENABLE_PYTHON "Enable Python bindings" ON)

option(SP2_BUILD_LIB     "Build library instead of binary" OFF)
option(BUILD_SHARED_LIBS  "Build shared libraries (DLLs)." OFF)

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
## Source files                                                               ##
################################################################################

set(SOURCE_FILES
        src/common/function_types.hpp
        src/common/function_types.cpp
        src/common/structure_t.cpp
        src/common/structure_t.hpp
        src/common/atom_types.hpp
        src/common/math/vec3_t.hpp
        src/common/math/vec3_util.hpp
        src/common/math/rotations.cpp
        src/common/math/rotations.hpp
        src/common/math/blas.hpp
        src/common/math/numerical_diff.hpp
        src/common/math/interpolate.cpp
        src/common/math/misc.hpp
        src/common/math/mat3x3_t.hpp
        src/common/util/random.cpp
        src/common/util/random.hpp
        src/common/util/templates.hpp
        src/common/util/timing.hpp
        src/common/util/timing.cpp
        src/common/util/modeling.cpp
        src/common/util/modeling.hpp
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
        src/common/io/gnuplot/gplot_structure.cpp
        src/common/io/gnuplot/gplot_structure.hpp
        src/common/json/json.cpp
        src/common/json/json.hpp
        src/common/json/json_serializable_t.hpp
        src/common/json/json_serializable_t.cpp
        src/common/minimize/minimize.hpp
        src/common/minimize/settings.cpp
        src/common/minimize/settings.hpp
        src/common/minimize/test/test_functions.cpp
        src/common/minimize/test/test_functions.hpp
        src/common/minimize/linesearch.cpp
        src/common/minimize/linear_cg.cpp
        src/common/minimize/acgsd.cpp
        src/common/minimize/metropolis.cpp
        src/common/minimize/fire.cpp
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
        IF_PHONOPY_ENABLED src/phonopy/bond_polarization.cpp
        IF_PHONOPY_ENABLED src/phonopy/bond_polarization.hpp
        IF_PHONOPY_ENABLED src/phonopy/phonopy_io.cpp
        IF_PHONOPY_ENABLED src/phonopy/phonopy_io.hpp
################################################################################
        IF_TESTS_ENABLED src/common/graph/test/test_ud_graph_t.cpp
        src/common/graph/ud_graph_t.cpp
        src/common/graph/ud_graph_t.hpp
        src/common/graph/permutation_t.cpp
        src/common/graph/permutation_t.hpp
        src/common/graph/graph.cpp
        src/common/graph/graph.hpp
        src/common/graph/directed_graph_t.cpp
        src/common/graph/directed_graph_t.hpp
        src/common/neighbor/utility_functions.cpp
        src/common/neighbor/utility_functions.hpp
        src/common/neighbor/periodic_cell_t.hpp
        src/common/neighbor/cell_array_t.cpp
        src/common/neighbor/cell_array_t.hpp
        src/common/neighbor/bond_control_t.cpp
        src/common/neighbor/bond_control_t.hpp
        src/common/enums.hpp
        src/common/neighbor/neighbor_control_t.hpp
        src/common/neighbor/neighbor_control_t.cpp
        src/common/env/env_t.cpp
        src/common/env/env_t.hpp
        src/airebo/system_control_t.cpp
        src/airebo/system_control_t.hpp
        src/airebo/interpolation_coeff.hpp
        src/airebo/utility_functions.hpp
        src/airebo/test.hpp
        src/airebo/test.cpp
        src/airebo/airebo_util.hpp
        # src/atac/system_control_t.cpp
        # src/atac/system_control_t.hpp
        src/atac/settings_t.hpp
        src/atac/settings_t.cpp
        src/run/common.cpp
        src/run/run_types.hpp
        src/run/run_settings_t.cpp
        src/run/run_settings_t.hpp
        src/run/program_args_t.cpp
        src/run/program_args_t.hpp
        src/run/run_relaxation.cpp
        # src/run/run_atac.cpp
        src/relaxation/relaxation_settings.cpp
        src/relaxation/relaxation_settings.hpp
        src/include/sp2/sp2.h src/common/python/bindings.cpp src/common/python/bindings.h src/common/util/unix.cpp src/common/util/unix.h)

# enable/disable inclusion of source files defined in ${SOURCE_FILES}
# depending on current options (e.g. SP2_ENABLE_MPI, SP2_ENABLE_PHONOPY, etc)
include(Sp2ParseSources)
parse_sources(SOURCE_FILES ${SOURCE_FILES})

################################################################################
## Main targets                                                               ##
################################################################################

if(NOT SP2_BUILD_LIB)
    # we are building an executable
    add_executable(sp2
            src/main.cpp
            ${SOURCE_FILES})
else()
    # we are building a library
    add_library(sp2
            src/run/library_interface.cpp
            ${SOURCE_FILES})

    target_include_directories(sp2 PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
            $<INSTALL_INTERFACE:include>)

    # install library headers
    install(DIRECTORY src/include/ DESTINATION include)
endif()

target_include_directories(sp2 PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR})

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    target_compile_definitions(sp2 PRIVATE
            SP2_DEBUG)
else()
    target_compile_options(sp2 PRIVATE
            -march=native)
endif()

install(TARGETS sp2
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

set_target_properties(sp2 PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF)

include(Sp2CompileChecks)

################################################################################
## Dependencies/external code                                                 ##
################################################################################

# add subdirectories for external libraries/code
add_subdirectory(external)

if(SP2_ENABLE_TESTS)
    add_test(all sp2 --test)
    target_link_libraries(sp2 PRIVATE gtest)
    # note: should work for gcc, clang
    target_compile_options(sp2 PRIVATE -Wall -Wextra -pedantic)
endif()

if(SP2_ENABLE_LAMMPS)
    target_link_libraries(sp2 PRIVATE lammps)
endif()

target_link_libraries(sp2 PRIVATE
        cxxopts
        jsoncpp_lib_static
        yaml-cpp
        templ
        randutils)

# internal code generation target
add_subdirectory(src/common/codegen)
add_dependencies(sp2 code_generation)

################################################################################
## MPI + boost                                                                ##
################################################################################

if(SP2_ENABLE_MPI)
    target_compile_definitions(sp2 PRIVATE
            SP2_ENABLE_MPI)

    # System MPI
    find_package(MPI REQUIRED)
    list(APPEND CMAKE_CXX_COMPILE_FLAGS
            ${MPI_CXX_COMPILE_FLAGS})

    set(CMAKE_CXX_LINK_FLAGS
            ${CMAKE_CXX_LINK_FLAGS}
            ${MPI_CXX_LINK_FLAGS})

    # Boost MPI (also requires Boost serialization)
    if(NOT BUILD_SHARED_LIBS)
        set(Boost_USE_STATIC_LIBS ON)
    endif()

    set(Boost_USE_MULTITHREADED ON)
    find_package(Boost REQUIRED COMPONENTS mpi serialization)

    # update target
    target_include_directories(sp2 PRIVATE
            ${MPI_CXX_INCLUDE_PATH}
            ${Boost_INCLUDE_DIRS})

    target_link_libraries(sp2 PUBLIC
            ${MPI_CXX_LIBRARIES}
            ${Boost_LIBRARIES})
else()
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)

    target_link_libraries(sp2
            Threads::Threads)
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

if (SP2_ENABLE_PYTHON)
    find_package(PythonLibs 3 REQUIRED)
    target_link_libraries(sp2 PRIVATE ${PYTHON_LIBRARIES})
    target_include_directories(sp2 PRIVATE ${PYTHON_INCLUDE_DIRS})
    target_compile_definitions(sp2 PRIVATE SP2_ENABLE_PYTHON)
endif()

if(SP2_ENABLE_PHONOPY)
    target_compile_definitions(sp2 PRIVATE
            SP2_ENABLE_PHONOPY)

    if(NOT PYMOD_PHONOPY_LOCATION)
        # first we try to find the python 3+ version of phonopy
        find_phonopy(3)
        # then if that fails, fall back to 2.7
        if(NOT PYMOD_PHONOPY_FOUND)
            find_phonopy(2.7 REQUIRED)
        endif()
    endif()
endif()


