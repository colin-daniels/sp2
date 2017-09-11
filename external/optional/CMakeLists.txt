download_project(PROJ optionaldl
    GIT_REPOSITORY https://github.com/akrzemi1/Optional
    GIT_TAG master
    UPDATE_DISCONNECTED 1
    QUIET)
add_library(optional INTERFACE)
target_include_directories(optional INTERFACE
    $<BUILD_INTERFACE:${optionaldl_SOURCE_DIR}>)
#download_project(PROJ optional
#    GIT_REPOSITORY https://github.com/akrzemi1/Optional
#    GIT_TAG master
#    UPDATE_DISCONNECTED 1
#    QUIET)

#target_link_libraries(mytarget ... akrzemi1::optional ...)
#add_library(optional INTERFACE)
#add_subdirectory(${optional_SOURCE_DIR} ${optional_BINARY_DIR})
#target_link_libraries(optional INTERFACE)

# ERROR: Missing an argument  (the BINARY_DIR is probably unset?)
#add_subdirectory(${optional_SOURCE_DIR})# ${optional_BINARY_DIR})

# Tried this but #include still can't find <optional.hpp>
#target_include_directories(optional INTERFACE
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
