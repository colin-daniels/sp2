include(DownloadProject)

# Note: if new directories are added in libs/, cmake needs to be re-run
# manually, as it will not know that the new directory is a dependency
file(GLOB LIB_SUBDIRS RELATIVE
        ${CMAKE_CURRENT_LIST_DIR}
        ${CMAKE_CURRENT_LIST_DIR}/*)

foreach(child ${LIB_SUBDIRS})
    if(IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
        add_subdirectory(${child})
    endif()
endforeach()
