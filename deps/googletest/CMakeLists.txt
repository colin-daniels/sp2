if(SP2_ENABLE_TESTS)
    # don't add gtest if a target already exists
    if(NOT TARGET gtest)
        # BSD 3-Clause: https://github.com/google/googletest/blob/release-1.10.0/LICENSE
        download_project(PROJ googletest
                GIT_REPOSITORY https://github.com/google/googletest.git
                GIT_TAG release-1.10.0
                UPDATE_DISCONNECTED 1
                QUIET)

        option(BUILD_GTEST "" ON)
        option(BUILD_GMOCK "" OFF)
        add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif()

    target_link_libraries(project-deps INTERFACE
            gtest)
endif()
