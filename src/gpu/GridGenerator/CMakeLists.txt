

vf_add_library(BUILDTYPE shared PRIVATE_LINK basics)

linkCUDA()
linkOpenMP()

vf_get_library_name(library_name)
set_target_properties(${library_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)