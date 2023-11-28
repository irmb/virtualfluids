PROJECT(SphereInChannel LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator FILES SphereInChannel.cpp)

set_source_files_properties(SphereInChannel.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(SphereInChannel PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
