

vf_add_library(BUILDTYPE shared DEPENDS basics)

include (${VF_CMAKE_DIR}/3rd/cuda.cmake)
include (${VF_CMAKE_DIR}/3rd/OpenMP.cmake)

vf_get_library_name(library_name)

set_property(TARGET ${library_name} PROPERTY CUDA_SEPARABLE_COMPILATION ON)