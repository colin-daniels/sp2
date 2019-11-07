cmake_minimum_required(VERSION 3.0)
project(randutils CXX)

add_library(randutils INTERFACE)

target_include_directories(randutils INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(project-deps INTERFACE
        randutils)
