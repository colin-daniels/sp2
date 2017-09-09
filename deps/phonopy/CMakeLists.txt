if(NOT SP2_ENABLE_PHONOPY OR PYMOD_PHONOPY_LOCATION)
    return()
endif()

include(FindPythonModule)

# force find_package() to look for python even if it already found it
unset(PYTHON_EXECUTABLE CACHE)
unset(PYTHONINTERP_FOUND CACHE)

# note: there's no library or target to be linked here, cmake will simply
# fail if it cannot find phonopy
find_package(PythonInterp 3 REQUIRED)
find_python_module(phonopy REQUIRED)
