PROJECT(GridGeneratorTest LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES gridGenerator.cpp)

set_source_files_properties(DrivenCavity.cpp PROPERTIES LANGUAGE CUDA)
