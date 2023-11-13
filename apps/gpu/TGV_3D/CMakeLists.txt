PROJECT(TGV_3D LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator FILES TGV_3D.cpp)

set_source_files_properties(TGV_3D.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(TGV_3D PROPERTIES CUDA_SEPARABLE_COMPILATION ON)