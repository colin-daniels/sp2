download_project(PROJ yamlcpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET)

set(YAML_CPP_BUILD_TOOLS OFF)
add_subdirectory(${yamlcpp_SOURCE_DIR} ${yamlcpp_BINARY_DIR})