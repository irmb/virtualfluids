


vf_get_library_name(library_name)
find_package(MPI REQUIRED)
target_include_directories(${library_name} PRIVATE ${MPI_C_INCLUDE_PATH})

target_link_libraries(${library_name} ${MPI_C_LIBRARIES})

