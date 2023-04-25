PROJECT(TGV_3D_GridRef LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator FILES TGV_3D_GridRef.cpp)

set_source_files_properties(TGV_3D_GridRef.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(TGV_3D_GridRef PROPERTIES CUDA_SEPARABLE_COMPILATION ON)