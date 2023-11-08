PROJECT(GridGeneratorTest LANGUAGES CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator MPI::MPI_CXX FILES gridGenerator.cpp)
