PROJECT(SphereMultiGPU LANGUAGES CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator FILES SphereMultiGPU.cpp)
