cmake_minimum_required(VERSION 3.2)

project(sp2 CXX)

# add extra modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
include(PreventInSourceBuilds)
include(DownloadProject)

# enable ctest
enable_testing()

option(SP2_USE_MPI "Enable MPI" ON)
option(SP2_USE_LAMMPS "Enable LAMMPS" OFF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set default build type to release
if(NOT WIN32)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: \n\
None Debug Release RelWithDebInfo MinSizeRel Coverage."
            FORCE)
    endif()
endif()

# add external code, sets SOURCE_FILES and EXTERNAL_INCLUDE_DIRS
add_subdirectory(external)
include_directories(src ${EXTERNAL_INCLUDE_DIRS})

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_definitions(-DSP2_DEBUG)

    add_test(all sp2 --test)

    link_libraries(gtest)
endif()

################################################################################
## Code Generation                                                            ##
################################################################################

# set output directory for generated code, make it, and add its parent
# directory to includes
set(CG_OUTPUT_DIR ${CMAKE_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${CG_OUTPUT_DIR})
include_directories(${CG_OUTPUT_DIR}/../)

# load the codegen project, sets CG_OUTPUT_FILES
add_subdirectory(src/common/codegen)

# custom command/target to run code generation
add_custom_command(OUTPUT ${CG_OUTPUT_FILES}
        COMMAND codegen
        DEPENDS codegen
        WORKING_DIRECTORY ${CG_OUTPUT_DIR}
        COMMENT "Generating source code")

add_custom_target(generate DEPENDS ${CG_OUTPUT_FILES}
        COMMENT "Verifying generated code")

################################################################################
## MPI                                                                        ##
################################################################################

if(SP2_USE_MPI)
    add_definitions(-DSP2_USE_MPI)

    find_package(MPI REQUIRED)
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

    link_libraries(SP2_LIBS ${SP2_LIBS}
            Threads::Threads)
endif()

################################################################################
## Main targets                                                               ##
################################################################################

set(SOURCE_FILES ${SOURCE_FILES}
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
        src/common/util/mpi.hpp
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
        src/common/minimize/pso/adaptive_pso.cpp
        src/common/minimize/pso/particle_t.cpp
        src/common/minimize/pso/particle_t.hpp
        src/common/minimize/pso/swarm_t.cpp
        src/common/minimize/pso/swarm_t.hpp
        src/common/minimize/acgsd.cpp
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
        src/symm/pso_adapters.cpp
        src/symm/pso_adapters.hpp
        src/symm/util.cpp
        src/symm/util.hpp
        src/symm/space_group_t.cpp
        src/symm/space_group_t.hpp
        src/symm/symm_settings_t.cpp
        src/symm/symm_settings_t.hpp
        src/symm/system_control_t.cpp
        src/symm/system_control_t.hpp
        # src/atac/system_control_t.cpp
        # src/atac/system_control_t.hpp
        src/atac/settings_t.hpp
        src/atac/settings_t.cpp
        # src/lammps/system_control_t.cpp
        # src/lammps/lammps_interface.hpp
        src/lammps/settings_t.cpp
        src/lammps/settings_t.hpp
        src/run/common.cpp
        src/run/run_symm.cpp
        # src/run/run_relaxation.cpp
        # src/run/run_atac.cpp
        src/run/run_types.hpp
        src/common/io/file_types/poscar.cpp
        src/common/io/file_types/poscar.hpp
        src/run/run_phonopy.cpp
        src/run/run_settings_t.cpp
        src/run/run_settings_t.hpp
        src/phonopy/phonopy_settings.cpp
        src/phonopy/phonopy_settings.hpp
        src/common/json/json_serializable_t.cpp
        src/airebo/airebo_util.hpp)

# sp2 executable (sp2)
add_executable(sp2-bin
        src/main.cpp
        ${SOURCE_FILES})

# explicitly set executable filename to 'sp2'
set_target_properties(sp2-bin PROPERTIES
        OUTPUT_NAME sp2)

# sp2 library (libsp2.so)
add_library(sp2-lib SHARED
        src/run/library_interface.cpp
        ${SOURCE_FILES})

# explicitly set library filename to libsp2.so
# additionally, stop the library from being made by default
set_target_properties(sp2-lib PROPERTIES
        OUTPUT_NAME sp2
        EXCLUDE_FROM_ALL 1
        EXCLUDE_FROM_DEFAULT_BUILD 1)

# code generation dependencies
add_dependencies(sp2-bin generate)
add_dependencies(sp2-lib generate)

