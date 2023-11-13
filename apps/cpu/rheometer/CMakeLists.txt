########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(rheometer)

vf_add_library(BUILDTYPE binary PRIVATE_LINK cpu_core basics ${MPI_CXX_LIBRARIES} NonNewtonianFluids FILES rheometer.cpp )