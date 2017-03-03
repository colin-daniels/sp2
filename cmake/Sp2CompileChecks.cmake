include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS
        ${CMAKE_REQUIRED_FLAGS}
        ${CMAKE_CXX14_STANDARD_COMPILE_OPTION}
        ${CMAKE_CXX_FLAGS})

set(CMAKE_REQUIRED_INCLUDES
        ${CMAKE_REQUIRED_INCLUDES}
        ${CMAKE_CURRENT_SOURCE_DIR}/include)

check_cxx_source_compiles("
    #include <unordered_map>
    enum class e : int {};
    std::unordered_map<e, int> m;

    int main() {}
" SP2_CAN_HASH_ENUM)

if(SP2_CAN_HASH_ENUM)
    target_compile_definitions(sp2 PUBLIC
            SP2_CAN_HASH_ENUM)
endif()