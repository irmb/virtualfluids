PROJECT(LBMAPP LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK lbmStatic FILES main.cpp)

set_source_files_properties(main.cpp PROPERTIES LANGUAGE CUDA)


set_target_properties(lbmApp PROPERTIES CUDA_SEPARABLE_COMPILATION ON)