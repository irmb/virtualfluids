function (linkMPI)

    find_package(MPI REQUIRED)

    vf_get_library_name(library_name)
    #target_include_directories(${library_name} PUBLIC ${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(${library_name} PRIVATE MPI::MPI_CXX)

endfunction()