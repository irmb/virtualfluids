PROJECT(TGV_3D_MultiGPU)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator FILES TGV_3D_MultiGPU.cpp)

set_source_files_properties(TGV_3D_MultiGPU.cpp PROPERTIES LANGUAGE CUDA)

linkCUDA()
