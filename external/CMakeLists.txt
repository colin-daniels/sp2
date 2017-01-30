set(SOURCE_FILES ${SOURCE_FILES}
        ${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp/jsoncpp.cpp
        PARENT_SCOPE)

set(EXTERNAL_INCLUDE_DIRS ${EXTERNAL_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp)

add_subdirectory(templ)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_subdirectory(googletest)
endif()

if(SP2_USE_LAMMPS)
    add_subdirectory(lammps)
endif()

# propagate include directories up
set(EXTERNAL_INCLUDE_DIRS ${EXTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
