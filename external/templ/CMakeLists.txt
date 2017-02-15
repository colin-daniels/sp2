download_project(PROJ templ
    GIT_REPOSITORY https://github.com/colin-daniels/templ.git
    GIT_TAG master
    UPDATE_DISCONNECTED 1
    QUIET)

set(TEMPL_ENABLE_INSTALL OFF)
set(TEMPL_ENABLE_TESTS OFF)
add_subdirectory(${templ_SOURCE_DIR} ${templ_BINARY_DIR})