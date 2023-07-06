PROJECT(MusselOyster LANGUAGES CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES MusselOyster.cpp)
