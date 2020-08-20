
find_package(OpenMP)
if (OPENMP_FOUND)
	#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    #set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler ${OpenMP_CXX_FLAGS}")
endif()

#message (Cuda Flags: ${OpenMP_CXX_FLAGS})
vf_get_library_name(library_name)
if(OpenMP_CXX_FOUND)
	target_link_libraries(${library_name} PUBLIC OpenMP::OpenMP_CXX)
endif()
