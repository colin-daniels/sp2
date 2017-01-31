# MODULE:   FindPythonModule
#
# PROVIDES:
#   find_python_module(<module> [REQUIRED])
#
#===============================================================================
#
# This code requires PYTHON_EXEC to be set and sets the following variables:
#
# ::
#
#   PYMOD_<uppercase module>_FOUND         - Was the Python module found
#   PYMOD_<uppercase module>_LOCATION      - path to the Python module
#
#===============================================================================
# Credit to Mark Moll, mmol@rice.edu
# Source: http://public.kitware.com/pipermail/cmake/2011-January/041666.html
#
# Modified by Colin Daniels
#===============================================================================

include(CMakeParseArguments)

function(find_python_module module)
    set(options REQUIRED)
    set(oneValueArgs "")
    set(multiValueArgs "")

    cmake_parse_arguments(FPM_ARGS
            "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT PYTHON_EXECUTABLE)
        message(FATAL_ERROR "find_python_module called without "
                "PYTHON_EXECUTABLE set.")
    endif()

    # convert module argument to package name
    string(TOUPPER ${module} package_name)
    set(package_name PYMOD_${package_name})

    if(FPM_ARGS_REQUIRED)
        set(${package_name}_FIND_REQUIRED TRUE)
    endif()

    if(NOT ${package_name}_FOUND)
        # A module's location is usually a directory, but for binary modules
        # it's a .so file.
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
                "import re, ${module}; print(re.compile('/__init__.py.*').sub('',${module}.__file__))"
                RESULT_VARIABLE _${module}_status
                OUTPUT_VARIABLE _${module}_location
                ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT _${module}_status)
            set(${package_name}_LOCATION ${_${module}_location} CACHE STRING
                    "Location of Python module ${module}")
        endif()

        find_package_handle_standard_args(${package_name}
                DEFAULT_MSG ${package_name}_LOCATION)
    endif()
endfunction(find_python_module)