project(lbm LANGUAGES CXX)

vf_add_library(NAME lbm PUBLIC_LINK basics)

if(BUILD_VF_GPU)
    add_subdirectory(cuda)
endif()


vf_add_tests()

#set_target_properties(lbm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
