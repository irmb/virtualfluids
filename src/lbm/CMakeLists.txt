
vf_add_library(PUBLIC_LINK basics)

if(BUILD_VF_GPU OR BUILD_VF_GKS)
    set_target_properties(lbm PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)

    set_source_files_properties(KernelParameter.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(CumulantChimera.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(BGK.cpp PROPERTIES LANGUAGE CUDA)
endif()

vf_add_tests()