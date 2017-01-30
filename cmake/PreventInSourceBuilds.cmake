#
# This function will prevent in-source builds
function(AssureOutOfSourceBuilds)
  # make sure the user doesn't play dirty with symlinks
  get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

  # disallow in-source builds
  if("${srcdir}" STREQUAL "${bindir}")
    message("######################################################")
    message("# This project should not be configured & built in its source directory")
    message("# You must run cmake in a build directory, e.g. the existing directory 'build'.")
    message("#")
    message("# You can proceed to configure and build by using the following commands")
    message("#")
    message("# cd build")
    message("# cmake ..")
    message("# make")
    message("#")
    message("# NOTE: Given that you already tried to make an in-source build")
    message("#       CMake have already created several files & directories")
    message("#       in your source tree. Remove them by doing:")
    message("#")
    message("#       rm -r CMakeFiles CMakeCache.txt")
    message("#")
    message("######################################################")
    message(FATAL_ERROR "Quitting configuration")
  endif()
endfunction()

AssureOutOfSourceBuilds()
