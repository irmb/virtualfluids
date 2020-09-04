

vf_add_library(BUILDTYPE shared DEPENDS basics)

linkCUDA()
linkOpenMP()

set_property(TARGET ${library_name} PROPERTY CUDA_SEPARABLE_COMPILATION ON)