PROJECT(SphereRefined LANGUAGES CUDA CXX)

#LIST(APPEND CS_COMPILER_FLAGS_CXX "-DOMPI_SKIP_MPICXX" )

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES SphereRefined.cpp)

set_source_files_properties(SphereRefined.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(SphereRefined PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
