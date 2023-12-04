
vf_add_library(PUBLIC_LINK basics logger PRIVATE_LINK project_warnings project_options)

if(VF_ENABLE_GPU)
    set_target_properties(lbm PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)
endif()

vf_add_tests()
