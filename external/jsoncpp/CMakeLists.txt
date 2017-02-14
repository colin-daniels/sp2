#add_library(jsoncpp INTERFACE)
#target_include_directories(jsoncpp INTERFACE
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
#target_sources(jsoncpp INTERFACE
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp.cpp>)

download_project(PROJ jsoncpp
        GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET)

add_subdirectory(${jsoncpp_SOURCE_DIR} ${jsoncpp_BINARY_DIR})