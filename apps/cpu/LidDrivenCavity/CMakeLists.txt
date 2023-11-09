
PROJECT(LidDrivenCavity)

vf_add_library(NAME LidDrivenCavityCPU BUILDTYPE binary PRIVATE_LINK cpu_core FILES LidDrivenCavity.cpp)
