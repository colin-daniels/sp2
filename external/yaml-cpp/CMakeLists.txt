download_project(PROJ yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET)

set(YAML_CPP_BUILD_TOOLS OFF)
add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR}
        EXCLUDE_FROM_ALL)