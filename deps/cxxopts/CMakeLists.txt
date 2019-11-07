# MIT License: https://github.com/jarro2783/cxxopts/blob/v2.2.0/LICENSE
download_project(PROJ cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG v2.2.0
        UPDATE_DISCONNECTED 1
        QUIET)

option(CXXOPTS_BUILD_EXAMPLES "" OFF)
option(CXXOPTS_BUILD_TESTS "" OFF)
add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR})

target_link_libraries(project-deps INTERFACE
        cxxopts)
