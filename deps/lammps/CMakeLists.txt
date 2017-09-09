if(NOT SP2_ENABLE_LAMMPS)
    return()
endif()

set(LAMMPS_SOURCE_DIR ${CMAKE_BINARY_DIR}/lammps-build)
set(LAMMPS_BUILD_DIR ${LAMMPS_SOURCE_DIR}/lammps)

# make output dir
file(MAKE_DIRECTORY ${LAMMPS_SOURCE_DIR})

# switch lammps library build mode based on static/shared
if(BUILD_SHARED_LIBS)
    set(LAMMPS_BUILD_MODE omp mode=shlib)
    set(LAMMPS_LIB_NAME liblammps.so)
else()
    set(LAMMPS_BUILD_MODE omp mode=lib)
    set(LAMMPS_LIB_NAME liblammps.a)
endif()

find_package(OpenMP REQUIRED)
ExternalProject_Add(lammps_external
        URL http://lammps.sandia.gov/tars/lammps-stable.tar.gz
        SOURCE_DIR "${LAMMPS_SOURCE_DIR}"
        CONFIGURE_COMMAND ""
        BINARY_DIR "${LAMMPS_BUILD_DIR}"
        # we use $(MAKE) to propagate -jN
        BUILD_COMMAND $(MAKE) yes-user-omp &&
            $(MAKE) --silent ${LAMMPS_BUILD_MODE}
        INSTALL_COMMAND ""
        STEP_TARGETS pre_config)

# custom config step to patch lammps
ExternalProject_Add_Step(lammps_external pre_config
        COMMAND ${CMAKE_COMMAND} -E rename "${LAMMPS_SOURCE_DIR}/src"
            "${LAMMPS_BUILD_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/CH.airebo.hpp"
            "${LAMMPS_BUILD_DIR}/"
        COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/pair_airebo.cpp"
            "${LAMMPS_BUILD_DIR}/"
        COMMAND sed -i "s/\\(CCFLAG.\\+\\s\\)-restrict/\\1/"
            "${LAMMPS_BUILD_DIR}/MAKE/OPTIONS/Makefile.omp"
        DEPENDERS configure
        DEPENDEES patch
        DEPENDS CH.airebo.hpp pair_airebo.cpp
        WORKING_DIRECTORY "${LAMMPS_SOURCE_DIR}")

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
        INTERFACE_INCLUDE_DIRECTORIES "${LAMMPS_SOURCE_DIR}"
        INTERFACE_COMPILE_DEFINITIONS SP2_ENABLE_LAMMPS)

target_link_libraries(project-deps INTERFACE
        lammps)
