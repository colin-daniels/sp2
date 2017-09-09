################################################################################
## Source files                                                               ##
################################################################################

set(SOURCE_FILES
        common/function_types.hpp
        common/function_types.cpp
        common/structure_t.cpp
        common/structure_t.hpp
        common/atom_types.hpp
        common/math/vec3_t.hpp
        common/math/vec3_util.hpp
        common/math/rotations.cpp
        common/math/rotations.hpp
        common/math/blas.hpp
        common/math/numerical_diff.hpp
        common/math/interpolate.cpp
        common/math/misc.hpp
        common/math/mat3x3_t.hpp
        common/util/random.cpp
        common/util/random.hpp
        common/util/templates.hpp
        common/util/timing.hpp
        common/util/timing.cpp
        common/util/modeling.cpp
        common/util/modeling.hpp
        common/io/structure.cpp
        common/io/structure.hpp
        common/io/util.cpp
        common/io/util.hpp
        common/io/file_types/xyz.cpp
        common/io/file_types/xyz.hpp
        common/io/file_types/cif.cpp
        common/io/file_types/cif.hpp
        common/io/file_types/poscar.cpp
        common/io/file_types/poscar.hpp
        common/io/gnuplot/gplot_structure.cpp
        common/io/gnuplot/gplot_structure.hpp
        common/json/json.cpp
        common/json/json.hpp
        common/json/json_serializable_t.hpp
        common/json/json_serializable_t.cpp
        common/minimize/minimize.hpp
        common/minimize/settings.cpp
        common/minimize/settings.hpp
        common/minimize/test/test_functions.cpp
        common/minimize/test/test_functions.hpp
        common/minimize/linesearch.cpp
        common/minimize/linear_cg.cpp
        common/minimize/acgsd.cpp
        common/minimize/metropolis.cpp
        common/minimize/metropolis.hpp
        common/minimize/metropolis_enums.hpp
        common/minimize/fire.cpp
################################################################################
        diagnostic/expect_python_headers
        diagnostic/forward_dependency_cycle
################################################################################
        IF_PYTHON_ENABLED ___________________ common/python.hpp
        IF_PYTHON_ENABLED ___________________ common/python/internals.hpp

        IF_PYTHON_ENABLED ___________________ common/python/conversion.hpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/conversion/custom.cpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/conversion/numpy.cpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/conversion/primitive.cpp
        IF_PYTHON_ENABLED ___________________ common/python/conversion/base_monomorphic.hpp
        IF_PYTHON_ENABLED ___________________ common/python/conversion/base_monomorphic_fwd.hpp
        IF_PYTHON_ENABLED ___________________ common/python/conversion/base_generic.hpp
        IF_PYTHON_ENABLED ___________________ common/python/conversion/base_generic_fwd.hpp
        IF_PYTHON_ENABLED ___________________ common/python/conversion/base_generic_raw.hpp
        IF_PYTHON_ENABLED ___________________ common/python/conversion/base_generic_raw_fwd.hpp
        IF_PYTHON_ENABLED ___________________ common/python/conversion/concrete.hpp
        IF_PYTHON_ENABLED ___________________ common/python/conversion/numpy.hpp

        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/environment.cpp
        IF_PYTHON_ENABLED ___________________ common/python/environment.hpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/error.cpp
        IF_PYTHON_ENABLED ___________________ common/python/error.hpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/include_numpy.cpp
        IF_PYTHON_ENABLED ___________________ common/python/include_numpy.hpp

        IF_PYTHON_ENABLED ___________________ common/python/modules/ext_modules.hpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/modules/ext_modules.cpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/modules/ext_modules/example.cpp

        IF_PYTHON_ENABLED ___________________ common/python/modules/fake_modules.hpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/modules/fake_modules.cpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/modules/fake_modules/example.cpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/modules/fake_modules/mutation_helper.cpp

        IF_PYTHON_ENABLED ___________________ common/python/types/as_ndarray.hpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/types/py_ref_t.cpp
        IF_PYTHON_ENABLED ___________________ common/python/types/py_ref_t.hpp
        IF_PYTHON_ENABLED ___________________ common/python/types/py_ref_t_body.hpp
        IF_PYTHON_ENABLED ___________________ common/python/types/py_ref_t_fwd.hpp
        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/types/py_object_t.cpp
        IF_PYTHON_ENABLED ___________________ common/python/types/py_object_t.hpp
        IF_PYTHON_ENABLED ___________________ common/python/types/py_object_t_body.hpp
        IF_PYTHON_ENABLED ___________________ common/python/types/py_object_t_fwd.hpp

        IF_PYTHON_ENABLED WITH_PYTHON_O_FLAGS common/python/utility.cpp
        IF_PYTHON_ENABLED ___________________ common/python/utility.hpp
################################################################################
        IF_MPI_ENABLED common/util/mpi.hpp
        # PSO requires MPI
        IF_MPI_ENABLED common/minimize/pso/adaptive_pso.cpp
        IF_MPI_ENABLED common/minimize/pso/particle_t.cpp
        IF_MPI_ENABLED common/minimize/pso/particle_t.hpp
        IF_MPI_ENABLED common/minimize/pso/swarm_t.cpp
        IF_MPI_ENABLED common/minimize/pso/swarm_t.hpp
        # and by extension so does symmetry search
        IF_MPI_ENABLED run/run_symm.cpp
        IF_MPI_ENABLED symm/pso_adapters.cpp
        IF_MPI_ENABLED symm/pso_adapters.hpp
        IF_MPI_ENABLED symm/util.cpp
        IF_MPI_ENABLED symm/util.hpp
        IF_MPI_ENABLED symm/space_group_t.cpp
        IF_MPI_ENABLED symm/space_group_t.hpp
        IF_MPI_ENABLED symm/symm_settings_t.cpp
        IF_MPI_ENABLED symm/symm_settings_t.hpp
        IF_MPI_ENABLED symm/system_control_t.cpp
        IF_MPI_ENABLED symm/system_control_t.hpp
################################################################################
        IF_LAMMPS_ENABLED lammps/system_control_t.cpp
        IF_LAMMPS_ENABLED lammps/lammps_interface.hpp
        IF_LAMMPS_ENABLED lammps/settings_t.cpp
        IF_LAMMPS_ENABLED lammps/settings_t.hpp
################################################################################
        IF_PHONOPY_ENABLED run/run_phonopy.cpp
        IF_PHONOPY_ENABLED phonopy/phonopy_settings.cpp
        IF_PHONOPY_ENABLED phonopy/phonopy_settings.hpp
        IF_PHONOPY_ENABLED phonopy/bond_polarization.cpp
        IF_PHONOPY_ENABLED phonopy/bond_polarization.hpp
        IF_PHONOPY_ENABLED phonopy/phonopy_io.cpp
        IF_PHONOPY_ENABLED phonopy/phonopy_io.hpp
################################################################################
        IF_TESTS_ENABLED common/graph/test/test_ud_graph_t.cpp
        common/graph/ud_graph_t.cpp
        common/graph/ud_graph_t.hpp
        common/graph/permutation_t.cpp
        common/graph/permutation_t.hpp
        common/graph/graph.cpp
        common/graph/graph.hpp
        common/graph/directed_graph_t.cpp
        common/graph/directed_graph_t.hpp
        common/neighbor/utility_functions.cpp
        common/neighbor/utility_functions.hpp
        common/neighbor/periodic_cell_t.hpp
        common/neighbor/cell_array_t.cpp
        common/neighbor/cell_array_t.hpp
        common/neighbor/bond_control_t.cpp
        common/neighbor/bond_control_t.hpp
        common/enums.hpp
        common/neighbor/neighbor_control_t.hpp
        common/neighbor/neighbor_control_t.cpp
        common/env/env_t.cpp
        common/env/env_t.hpp
        airebo/system_control_t.cpp
        airebo/system_control_t.hpp
        airebo/interpolation_coeff.hpp
        airebo/utility_functions.hpp
        airebo/test.hpp
        airebo/test.cpp
        airebo/airebo_util.hpp
        # atac/system_control_t.cpp
        # atac/system_control_t.hpp
        atac/settings_t.hpp
        atac/settings_t.cpp
        run/common.cpp
        run/run_types.hpp
        run/run_settings_t.cpp
        run/run_settings_t.hpp
        run/program_args_t.cpp
        run/program_args_t.hpp
        run/run_relaxation.cpp
        # run/run_atac.cpp
        relaxation/relaxation_settings.cpp
        relaxation/relaxation_settings.hpp
        include/sp2/sp2.h)

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
            main.cpp
            ${SOURCE_FILES})
else()
    # we are building a library
    add_library(sp2
            run/library_interface.cpp
            ${SOURCE_FILES})

    target_include_directories(sp2 PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)
endif()

set_target_properties(sp2 PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF)

include(Sp2CompileChecks)

# Note: project files should be included via something like:
# #include <src/subdir/file.h> where src is the name of the current list dir
get_filename_component(PARENT_DIR "${CMAKE_CURRENT_LIST_DIR}" DIRECTORY)
target_include_directories(sp2 PRIVATE BEFORE
        $<BUILD_INTERFACE:${PARENT_DIR}>
        # for generated headers
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
        # TODO: remove when done refactoring includes to use src/ prefix
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    target_compile_definitions(sp2 PRIVATE
            SP2_DEBUG)
else()
    target_compile_options(sp2 PRIVATE
            -march=native)
endif()

################################################################################
## Dependencies/external code                                                 ##
################################################################################

target_link_libraries(sp2 PRIVATE
        project-deps)

# internal code generation target
add_subdirectory(common/codegen)
add_dependencies(sp2 code_generation)

if(SP2_ENABLE_TESTS)
    add_test(all sp2 --test)
    target_compile_definitions(sp2 PRIVATE
            SP2_ENABLE_TESTS)

    # note: should work for gcc, clang
    target_compile_options(sp2 PRIVATE
            -Wall -Wextra -pedantic)
endif()

if(SP2_ENABLE_MPI)
    target_compile_definitions(sp2 PRIVATE
            SP2_ENABLE_MPI)
endif()

if (SP2_ENABLE_PYTHON)
    target_compile_definitions(sp2 PRIVATE
            SP2_ENABLE_PYTHON)
endif()

if(SP2_ENABLE_PHONOPY)
    target_compile_definitions(sp2 PRIVATE
            SP2_ENABLE_PHONOPY)
endif()

################################################################################
# Installation code                                                            #
################################################################################
# adapted mostly from github.com/forexample/package-example

# mostly install-related variables
set(config_install_dir "lib/cmake/${PROJECT_NAME}")
set(include_install_dir "include")
set(generated_dir "${PROJECT_BINARY_DIR}/generated")

# configuration
set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
set(targets_export_name "${PROJECT_NAME}Targets")
set(namespace "${PROJECT_NAME}::")

if(SP2_BUILD_LIB)
    # install library headers
    install(DIRECTORY include/ DESTINATION include)
endif()

# Include module with fuction 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
# Note: PROJECT_VERSION is used as a VERSION
message(STATUS "VERSION: ${PROJECT_NAME}_VERSION")
message(STATUS "VERSION: ${${PROJECT_NAME}_VERSION}")
write_basic_package_version_file("${version_config}"
        VERSION ${${PROJECT_NAME}_VERSION}
        COMPATIBILITY SameMajorVersion)

# Configure '<PROJECT-NAME>Config.cmake'
# Use variables:
#   * targets_export_name
#   * PROJECT_NAME
configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/Config.cmake.in"
        "${project_config}"
        INSTALL_DESTINATION "${config_install_dir}")

install(TARGETS sp2
        EXPORT ${targets_export_name}
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib"
        INCLUDES DESTINATION "${include_install_dir}")

# copy configs
install(FILES "${project_config}" "${version_config}"
        DESTINATION "${config_install_dir}")

install(EXPORT "${targets_export_name}"
        NAMESPACE "${namespace}"
        DESTINATION "${config_install_dir}")
