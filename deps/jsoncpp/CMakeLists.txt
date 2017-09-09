download_project(PROJ jsoncpp
        GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET)

option(JSONCPP_WITH_TESTS "" OFF)
option(JSONCPP_WITH_POST_BUILD_UNITTEST "" OFF)
option(JSONCPP_WITH_PKGCONFIG_SUPPORT "" OFF)
option(BUILD_SHARED_LIBS "" OFF)
option(BUILD_STATIC_LIBS "" ON)
add_subdirectory(${jsoncpp_SOURCE_DIR} ${jsoncpp_BINARY_DIR}
        EXCLUDE_FROM_ALL)

target_link_libraries(project-deps INTERFACE
        jsoncpp_lib_static)
