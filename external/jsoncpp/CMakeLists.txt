add_library(jsoncpp INTERFACE)
target_include_directories(jsoncpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_sources(jsoncpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp.cpp>)