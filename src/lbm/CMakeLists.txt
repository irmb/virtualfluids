project(lbm LANGUAGES CXX CUDA)

vf_add_library(NAME lbm PUBLIC_LINK basics)

vf_add_tests()

set_target_properties(lbm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
