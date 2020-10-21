

vf_add_library(BUILDTYPE shared PRIVATE_LINK basics)
vf_get_library_name(library_name)

linkCUDA()
linkOpenMP(${library_name})

set_target_properties(${library_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)