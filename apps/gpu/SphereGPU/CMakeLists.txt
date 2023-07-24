PROJECT(SphereGPU LANGUAGES CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator FILES Sphere.cpp)
