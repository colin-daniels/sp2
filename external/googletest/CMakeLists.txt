download_project(PROJ googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET)

option(BUILD_GTEST "" ON)
option(BUILD_GMOCK "" OFF)
add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
