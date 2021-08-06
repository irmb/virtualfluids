########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(viskomat2)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} FILES viskomat2.cpp )