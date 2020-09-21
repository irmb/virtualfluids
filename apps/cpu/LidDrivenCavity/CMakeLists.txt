
PROJECT(LidDrivenCavity)

vf_add_library(NAME LidDrivenCavityCPU BUILDTYPE binary DEPENDS VirtualFluidsCore basics muparser FILES LidDrivenCavity.cpp)
