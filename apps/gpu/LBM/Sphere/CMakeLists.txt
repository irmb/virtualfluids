PROJECT(DrivenCavity LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES Sphere.cpp)

set_source_files_properties(Sphere.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(Sphere PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
