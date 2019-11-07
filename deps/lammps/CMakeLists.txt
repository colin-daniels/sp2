if(NOT SP2_ENABLE_LAMMPS)
    return()
endif()

set(LAMMPS_SOURCE_DIR ${CMAKE_BINARY_DIR}/lammps-build)
set(LAMMPS_BUILD_DIR ${LAMMPS_SOURCE_DIR}/src)

# make output dir
file(MAKE_DIRECTORY ${LAMMPS_SOURCE_DIR})

# Convert cmake lists to things that can be passed to the compiler on the
# command line.
string(REPLACE ";" " " MPI_CXX_LIBRARIES_STR "${MPI_CXX_LIBRARIES}")
if(MPI_CXX_INCLUDE_PATH)
    # note: needs -I at the beginning since cmake lists dont start with ;
    string(REPLACE ";" " -I"
            MPI_CXX_INCLUDE_PATH_STR "-I${MPI_CXX_INCLUDE_PATH}")
endif()

# configure LAMMPS makefile patch
configure_file(makefile.patch.in makefile.patch @ONLY)

# switch lammps library build mode based on static/shared
if(BUILD_SHARED_LIBS)
    set(LAMMPS_BUILD_MODE omp mode=shlib)
    set(LAMMPS_LIB_NAME liblammps.so)
else()
    set(LAMMPS_BUILD_MODE omp mode=lib)
    set(LAMMPS_LIB_NAME liblammps.a)
endif()

find_package(OpenMP REQUIRED)
# GPL v2.0: https://github.com/lammps/lammps/blob/patch_22Sep2017/LICENSE
ExternalProject_Add(lammps_external
        URL https://github.com/lammps/lammps/archive/patch_22Sep2017.tar.gz
        SOURCE_DIR "${LAMMPS_SOURCE_DIR}"

        # patch LAMMPS makefile to pass along correct compiler and MPI
        # compilation flags
        PATCH_COMMAND patch "${LAMMPS_BUILD_DIR}/MAKE/OPTIONS/Makefile.omp"
            "${CMAKE_CURRENT_BINARY_DIR}/makefile.patch"
        # patch pair_airebo to read fake parameter "file"
        COMMAND patch "${LAMMPS_BUILD_DIR}/MANYBODY/pair_airebo.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/airebo_parameters.patch"
        # copy header that contains fake airebo parameter "file"
        COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/CH.airebo.hpp"
            "${LAMMPS_BUILD_DIR}"

        # enable openmp user package
        CONFIGURE_COMMAND $(MAKE) yes-user-omp

        BINARY_DIR "${LAMMPS_BUILD_DIR}"
        # we use $(MAKE) to propagate -jN
        BUILD_COMMAND $(MAKE) --silent ${LAMMPS_BUILD_MODE}

        # last step adds a symlink to allow us to include the lammps headers
        # via #include <lammps/X> instead of the conflicting #include <src/X>
        INSTALL_COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${LAMMPS_BUILD_DIR}"
            "${LAMMPS_SOURCE_DIR}/lammps")

# create imported library representing lammps, its location, and headers
if(BUILD_SHARED_LIBS)
    add_library(lammps SHARED IMPORTED GLOBAL)
else()
    add_library(lammps STATIC IMPORTED GLOBAL)
endif()

target_link_libraries(lammps INTERFACE "-fopenmp")
add_dependencies(lammps lammps_external)
set_target_properties(lammps PROPERTIES
        IMPORTED_LOCATION "${LAMMPS_BUILD_DIR}/${LAMMPS_LIB_NAME}"
        INTERFACE_INCLUDE_DIRECTORIES "${LAMMPS_SOURCE_DIR}")

target_link_libraries(project-deps INTERFACE
        lammps)
