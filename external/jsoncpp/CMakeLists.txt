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

## python needed for amagalmate.py
#if (NOT PYTHONINTERP_FOUND)
#    find_package(PythonInterp 3.4)
#    find_package(PythonInterp 2.6 REQUIRED)
#endif()
#
#set(jsoncpp_DIR ${jsoncpp_SOURCE_DIR}/dist)
#add_custom_command(OUTPUT ${jsoncpp_DIR}/jsoncpp.cpp
#        ${jsoncpp_DIR}/json/json.h
#        ${jsoncpp_DIR}/json/json-forwards.h
#        COMMAND ${PYTHON_EXECUTABLE} amalgamate.py
#        WORKING_DIRECTORY ${jsoncpp_SOURCE_DIR})
#
#add_library(jsoncpp INTERFACE)
#target_include_directories(jsoncpp INTERFACE
#        $<BUILD_INTERFACE:${jsoncpp_DIR}>)
#target_sources(jsoncpp INTERFACE
#        $<BUILD_INTERFACE:${jsoncpp_DIR}/jsoncpp.cpp>)