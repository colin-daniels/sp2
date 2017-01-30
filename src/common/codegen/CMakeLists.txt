cmake_minimum_required(VERSION 3.2)

project(codegen CXX)

# code generation source files
set(SOURCE_FILES
        main.cpp
        util/lookup_table.hpp
        util/lookup_table.cpp
        util/util.cpp
        util/util.hpp
        util/stream_io.cpp
        util/stream_io.hpp
        outputs/outputs.hpp
        outputs/exp_table.cpp
        ../util/random.cpp
        ../util/random.hpp
        ../util/timing.hpp
        ../util/timing.cpp)

# code generation output directory and files
set(OUTPUT_FILES
        exp_table.hpp)

add_executable(codegen ${SOURCE_FILES})
target_include_directories(codegen PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(codegen PRIVATE -O3)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_test(all codegen --test)
endif()

# prepend output directory to files
foreach(OUTPUT_FILE ${OUTPUT_FILES})
    set(MODIFIED ${MODIFIED} ${CG_OUTPUT_DIR}/${OUTPUT_FILE})
endforeach()

# propagate up
set(CG_OUTPUT_FILES ${MODIFIED} PARENT_SCOPE)