########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(JetBreakup)

vf_add_library(BUILDTYPE binary PRIVATE_LINK cpu_core basics ${MPI_CXX_LIBRARIES} MultiphaseFlow NonNewtonianFluids FILES JetBreakup.cpp )