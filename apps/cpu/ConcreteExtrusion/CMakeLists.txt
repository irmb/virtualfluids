########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(ConcreteExtrusion)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} MultiphaseFlow NonNewtonianFluids FILES ConcreteExtrusion.cpp )