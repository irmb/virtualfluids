########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(sphere LANGUAGES CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} FILES sphere.cpp )