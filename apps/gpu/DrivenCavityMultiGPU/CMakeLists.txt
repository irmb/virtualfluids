PROJECT(DrivenCavityMultiGPU LANGUAGES CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES DrivenCavityMultiGPU.cpp)
