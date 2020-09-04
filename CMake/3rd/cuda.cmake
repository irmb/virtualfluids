
function(linkCUDA)

    find_path(CUDA_CUT_INCLUDE_DIR
      helper_cuda.h
      PATHS "$ENV{NVCUDASAMPLES_ROOT}" "${NVCUDASAMPLES_ROOT}"
      PATH_SUFFIXES "common/inc"
      DOC "Location of helper_cuda.h"
      NO_DEFAULT_PATH
    )

    vf_get_library_name(library_name)
    target_include_directories(${library_name} PRIVATE ${CUDA_INCLUDE_DIRS} ${CUDA_CUT_INCLUDE_DIR})

    # set the following properties only for specific targets
    # set_property(TARGET ${targetName} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    # set_property(TARGET ${targetName} PROPERTY CUDA_64_BIT_DEVICE_CODE ON)
endfunction()