PROJECT(DrivenCavityMultiGPU LANGUAGES CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator FILES DrivenCavityMultiGPU.cpp)
