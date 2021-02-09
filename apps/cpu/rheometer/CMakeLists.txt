########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(rheometer)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} FILES rheometer.cpp )