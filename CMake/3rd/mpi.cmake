function (linkMPI)

    find_package(MPI REQUIRED)

    vf_get_library_name(library_name)
    target_link_libraries(${library_name} PUBLIC MPI::MPI_CXX)

endfunction()