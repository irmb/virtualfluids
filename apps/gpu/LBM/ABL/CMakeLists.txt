PROJECT(ABL LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES ABL.cpp)

set_source_files_properties(ABL.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(ABL PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
