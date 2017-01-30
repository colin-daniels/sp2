set(LAMMPS_SOURCE_DIR ${CMAKE_BINARY_DIR}/lammps)
set(LAMMPS_BUILD_DIR ${LAMMPS_SOURCE_DIR}/lammps)

# switch lammps library build mode based on static/shared
if(SP2_STATIC)
    set(LAMMPS_BUILD_MODE mpi mode=lib)
    set(LAMMPS_LIB_NAME liblammps_mpi.a)
    set(LAMMPS_LIB_TYPE STATIC)
else()
    set(LAMMPS_BUILD_MODE mpi mode=shlib)
    set(LAMMPS_LIB_NAME liblammps_mpi.so)
    set(LAMMPS_LIB_TYPE SHARED)
endif()

file(MAKE_DIRECTORY ${LAMMPS_SOURCE_DIR})
ExternalProject_Add(lammps_external
        URL ${CMAKE_CURRENT_SOURCE_DIR}/lammps.tar.gz
        URL_HASH SHA1=4d9ad513c1ed9eef78a524f4ee7408b75c610b0f
        SOURCE_DIR ${LAMMPS_SOURCE_DIR}
        CONFIGURE_COMMAND ""
        BINARY_DIR ${LAMMPS_BUILD_DIR}
        BUILD_COMMAND make -j4 ${LAMMPS_BUILD_MODE}
        INSTALL_COMMAND ""
        STEP_TARGETS pre_config)

# custom config step to patch lammps
ExternalProject_Add_Step(lammps_external pre_config
        COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/CH.airebo.hpp"
            "${LAMMPS_BUILD_DIR}/"
        COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/pair_airebo.cpp"
            "${LAMMPS_BUILD_DIR}/"
        DEPENDERS configure
        DEPENDEES patch
        DEPENDS CH.airebo.hpp pair_airebo.cpp)

# create imported library representing lammps, its location, and headers
add_library(lammps ${LAMMPS_LIB_TYPE} IMPORTED GLOBAL)
add_dependencies(lammps lammps_external)
set_target_properties(lammps PROPERTIES
        IMPORTED_LOCATION "${LAMMPS_BUILD_DIR}/${LAMMPS_LIB_NAME}"
        INTERFACE_INCLUDE_DIRECTORIES "${LAMMPS_SOURCE_DIR}")
