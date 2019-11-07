# MIT License: https://github.com/jbeder/yaml-cpp/blob/yaml-cpp-0.6.3/LICENSE
download_project(PROJ yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG yaml-cpp-0.6.3
        UPDATE_DISCONNECTED 1
        QUIET)

option(YAML_CPP_BUILD_TOOLS "" OFF)
add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR}
        EXCLUDE_FROM_ALL)

target_include_directories(yaml-cpp INTERFACE
        "$<BUILD_INTERFACE:${yaml-cpp_SOURCE_DIR}/include>")

target_link_libraries(project-deps INTERFACE
        yaml-cpp)
