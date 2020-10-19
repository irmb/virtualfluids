PROJECT(TGV_3D)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator FILES TGV_3D.cpp)

linkCUDA()
