function (linkOpenMP)

	find_package(OpenMP REQUIRED)

	if(OpenMP_CXX_FOUND)
		vf_get_library_name(library_name)
		target_link_libraries(${library_name} PUBLIC OpenMP::OpenMP_CXX)
	endif()

endfunction()