
PROJECT(LidDrivenCavity)

vf_add_library(NAME LidDrivenCavityCPU BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore FILES LidDrivenCavity.cpp)
