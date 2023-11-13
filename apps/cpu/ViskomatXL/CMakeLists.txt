########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(viskomat)

vf_add_library(BUILDTYPE binary PRIVATE_LINK cpu_core basics ${MPI_CXX_LIBRARIES} NonNewtonianFluids FILES viskomat.cpp )