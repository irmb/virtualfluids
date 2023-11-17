
vf_add_library(PUBLIC_LINK logger)

if(MSVC)
    target_link_libraries(parallel PRIVATE ws2_32)
endif()

# TODO: Make VirtualFluids work without MPI
# if(BUILD_USE_MPI)
    target_link_libraries(parallel PUBLIC MPI::MPI_CXX)
# endif()
