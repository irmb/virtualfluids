

PROJECT(LidDrivenCavity)


vf_add_library(BUILDTYPE binary PRIVATE_LINK basics GridGenerator VirtualFluids_GPU GksMeshAdapter GksGpu FILES LidDrivenCavity.cpp)
