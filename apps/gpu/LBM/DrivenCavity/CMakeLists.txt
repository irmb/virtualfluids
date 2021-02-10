PROJECT(DrivenCavity LANGUAGES CUDA CXX)

#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-DOMPI_SKIP_MPICXX" )

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES DrivenCavity.cpp)

set_source_files_properties(DrivenCavity.cpp PROPERTIES LANGUAGE CUDA)
