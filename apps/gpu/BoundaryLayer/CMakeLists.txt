PROJECT(BoundaryLayer LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES BoundaryLayer.cpp)

set_source_files_properties(BoundaryLayer.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(BoundaryLayer PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
