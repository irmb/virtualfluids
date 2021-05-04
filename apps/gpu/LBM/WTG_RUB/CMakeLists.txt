PROJECT(WTG_RUB LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES WTG_RUB.cpp)

set_source_files_properties(WTG_RUB.cpp PROPERTIES LANGUAGE CUDA)
