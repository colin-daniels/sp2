if (NOT SP2_ENABLE_PYTHON)
    return()
endif()

find_package(PythonLibs 3 REQUIRED)
find_package(NumPy REQUIRED)

target_link_libraries(project-deps INTERFACE
        ${PYTHON_LIBRARIES})

target_include_directories(project-deps INTERFACE
        ${PYTHON_INCLUDE_DIRS}
        ${PYTHON_NUMPY_INCLUDE_DIR})
