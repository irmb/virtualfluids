
vf_add_library(PUBLIC_LINK logger MPI::MPI_CXX basics)

if(MSVC)
    target_link_libraries(parallel PRIVATE ws2_32)
endif()
