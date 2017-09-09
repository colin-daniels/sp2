if(SP2_ENABLE_MPI)
    # System MPI
    find_package(MPI REQUIRED)
    list(APPEND CMAKE_CXX_COMPILE_FLAGS
            ${MPI_CXX_COMPILE_FLAGS})

    set(CMAKE_CXX_LINK_FLAGS
            ${CMAKE_CXX_LINK_FLAGS}
            ${MPI_CXX_LINK_FLAGS})

    # Boost MPI (also requires Boost serialization)
    if(NOT BUILD_SHARED_LIBS)
        set(Boost_USE_STATIC_LIBS ON)
    endif()

    set(Boost_USE_MULTITHREADED ON)
    find_package(Boost REQUIRED COMPONENTS mpi serialization)

    # update target
    target_include_directories(project-deps INTERFACE
            ${MPI_CXX_INCLUDE_PATH}
            ${Boost_INCLUDE_DIRS})

    target_link_libraries(project-deps INTERFACE
            ${MPI_CXX_LIBRARIES}
            ${Boost_LIBRARIES})
else()
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)

    target_link_libraries(project-deps INTERFACE
            Threads::Threads)
endif()
