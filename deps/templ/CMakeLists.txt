download_project(PROJ templ
    GIT_REPOSITORY https://github.com/colin-daniels/templ.git
    GIT_TAG 5c2dd1d36b107cbeb5fc46a99b5775bbee355cad
    UPDATE_DISCONNECTED 1
    QUIET)

option(TEMPL_ENABLE_INSTALL "" OFF)
option(TEMPL_ENABLE_TESTS "" OFF)
add_subdirectory(${templ_SOURCE_DIR} ${templ_BINARY_DIR})

target_link_libraries(project-deps INTERFACE
        templ)
