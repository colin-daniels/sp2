download_project(PROJ googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET)

set(BUILD_GTEST ON)
set(BUILD_GMOCK ON)
add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
