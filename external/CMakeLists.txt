if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_subdirectory(googletest)
endif()

add_subdirectory(jsoncpp)
add_subdirectory(lammps)
add_subdirectory(templ)