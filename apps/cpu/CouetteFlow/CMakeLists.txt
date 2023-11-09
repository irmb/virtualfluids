########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(CouetteFlow)

vf_add_library(BUILDTYPE binary PRIVATE_LINK cpu_core basics ${MPI_CXX_LIBRARIES} NonNewtonianFluids FILES cflow.cpp )