

vf_add_library(PRIVATE_LINK basics)
vf_get_library_name(library_name)

linkCUDA()
linkOpenMP(${library_name})


set_target_properties(${library_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# according to linker error when building static libraries.
# https://stackoverflow.com/questions/50033435/cmake-cuda-separate-compilation-static-lib-link-error-on-windows-but-not-on-ubun
if (NOT BUILD_SHARED_LIBRARY)
    set_target_properties(${library_name} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()