PROJECT(TGV_3D)

vf_add_library(BUILDTYPE binary DEPENDS basics VirtualFluids_GPU GridGenerator FILES TGV_3D.cpp)

linkCUDA()
