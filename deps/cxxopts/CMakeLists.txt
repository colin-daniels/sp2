download_project(PROJ cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET)

option(CXXOPTS_BUILD_EXAMPLES "" OFF)
option(CXXOPTS_BUILD_TESTS "" OFF)
add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR})

target_link_libraries(project-deps INTERFACE
        cxxopts)
