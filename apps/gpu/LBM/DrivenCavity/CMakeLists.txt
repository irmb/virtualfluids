PROJECT(DrivenCavity LANGUAGES CUDA CXX)

#LIST(APPEND CS_COMPILER_FLAGS_CXX "-DOMPI_SKIP_MPICXX" )

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES DrivenCavity.cpp)

set_source_files_properties(DrivenCavity.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(DrivenCavity PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

#target_compile_options(DrivenCavity PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-dlto>)

#target_link_options(DrivenCavity PRIVATE $<DEVICE_LINK:-dlto>)