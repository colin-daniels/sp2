if(SP2_ENABLE_TESTS)
    add_subdirectory(googletest)
    target_compile_definitions(gtest INTERFACE
            SP2_ENABLE_TESTS)
endif()

if(SP2_ENABLE_LAMMPS)
    add_subdirectory(lammps)
endif()

add_subdirectory(templ)
add_subdirectory(cxxopts)
add_subdirectory(jsoncpp)
add_subdirectory(yaml-cpp)
add_subdirectory(randutils)
