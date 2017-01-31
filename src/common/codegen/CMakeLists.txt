# code generation target
add_executable(code_generation
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

set_target_properties(code_generation PROPERTIES
        OUTPUT_NAME codegen)

target_include_directories(code_generation PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})

# requires templ
target_link_libraries(code_generation templ)

if(SP2_ENABLE_TESTS)
    add_test(all codegen --test)
    target_compile_definitions(code_generation PRIVATE
            SP2_ENABLE_TESTS)
    target_link_libraries(code_generation gtest)
endif()

# code generation root directory (generated files go in dir/generated/)
set(OUTPUT_ROOT "${CMAKE_BINARY_DIR}")
file(MAKE_DIRECTORY "${OUTPUT_ROOT}/generated")

# custom command/target to run code generation
add_custom_command(TARGET code_generation POST_BUILD
        COMMAND code_generation
        DEPENDS code_generation
        WORKING_DIRECTORY "${OUTPUT_ROOT}/generated"
        COMMENT "Generating source code")


