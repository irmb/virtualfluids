macro(linkCuda)
#input: -

    find_path(CUDA_CUT_INCLUDE_DIR
      helper_cuda.h
      PATHS "$ENV{NVCUDASAMPLES8_0_ROOT}"
      PATH_SUFFIXES "common/inc" "v8.0/common/inc"
      DOC "Location of helper_cuda.h"
      NO_DEFAULT_PATH
    )
    
    INCLUDE_DIRECTORIES( ${CUDA_INCLUDE_DIRS} ${CUDA_CUT_INCLUDE_DIR})

	set_property(TARGET ${targetName} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
	#set_property(TARGET ${targetName} PROPERTY CUDA_64_BIT_DEVICE_CODE ON)
	
    #set(CUDA_NVCC_FLAGS "-G" CACHE TYPE INTERNAL FORCE)
    #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -gencode arch=compute_13,code=sm_13 )
	
endmacro(linkCuda)
