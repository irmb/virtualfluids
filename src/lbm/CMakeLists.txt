project(lbm LANGUAGES CXX)

if(BUILD_VF_GPU)
    enable_language(CUDA)
endif()

vf_add_library(PUBLIC_LINK basics)

if(BUILD_VF_GPU)
set_target_properties(lbm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
