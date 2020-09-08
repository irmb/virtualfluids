PROJECT(TGV_3D_MultiGPU)

vf_add_library(BUILDTYPE binary DEPENDS basics VirtualFluids_GPU GridGenerator FILES TGV_3D_MultiGPU.cpp)

linkCUDA()
