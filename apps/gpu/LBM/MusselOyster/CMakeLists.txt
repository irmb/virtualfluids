PROJECT(DrivenCavity LANGUAGES CUDA CXX)

#LIST(APPEND CS_COMPILER_FLAGS_CXX "-DOMPI_SKIP_MPICXX" )

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES MusselOyster.cpp)

set_source_files_properties(MusselOyster.cpp PROPERTIES LANGUAGE CUDA)
