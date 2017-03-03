include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS
        ${CMAKE_REQUIRED_FLAGS}
        ${CMAKE_CXX14_STANDARD_COMPILE_OPTION}
        ${CMAKE_CXX_FLAGS})

set(CMAKE_REQUIRED_INCLUDES
        ${CMAKE_REQUIRED_INCLUDES}
        ${PROJECT_SOURCE_DIR}/sp2)

check_cxx_source_compiles("
    #include <unordered_map>
    enum class E : int {};
    std::unordered_map<E, int> m;

    int main() {}
" SP2_CAN_HASH_ENUM)

if(SP2_CAN_HASH_ENUM)
    target_compile_definitions(sp2 PUBLIC
            SP2_CAN_HASH_ENUM)
endif()