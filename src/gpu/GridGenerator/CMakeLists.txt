

vf_add_library(BUILDTYPE shared DEPENDS basics)

linkCUDA()
linkOpenMP()

vf_get_library_name(library_name)
set_target_properties(${library_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)