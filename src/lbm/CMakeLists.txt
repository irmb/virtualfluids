
vf_add_library(PUBLIC_LINK basics)

if(BUILD_VF_GPU)
    set_target_properties(lbm PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)
endif()

vf_add_tests()