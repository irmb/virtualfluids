

PROJECT(LidDrivenCavity)


vf_add_library(BUILDTYPE binary DEPENDS basics GridGenerator VirtualFluids_GPU GksMeshAdapter GksGpu FILES LidDrivenCavity.cpp)
