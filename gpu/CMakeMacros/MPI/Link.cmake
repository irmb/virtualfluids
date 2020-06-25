macro(linkMPI targetName)

	find_package(MPI REQUIRED)
	include_directories(${MPI_C_INCLUDE_PATH})

	target_link_libraries(${targetName} ${MPI_C_LIBRARIES})

endmacro(linkMPI)
