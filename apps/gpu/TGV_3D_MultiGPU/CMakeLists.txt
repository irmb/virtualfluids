PROJECT(TGV_3D_MultiGPU LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator FILES TGV_3D_MultiGPU.cpp)

set_source_files_properties(TGV_3D_MultiGPU.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(TGV_3D_MultiGPU PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
