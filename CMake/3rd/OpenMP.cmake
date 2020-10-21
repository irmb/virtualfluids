function (linkOpenMP targetName)

	if(NOT USE_OPENMP)
		return()
	endif()

	find_package(OpenMP REQUIRED)

	if(OpenMP_CXX_FOUND)
		target_link_libraries(${targetName} PUBLIC OpenMP::OpenMP_CXX)
	endif()

endfunction()