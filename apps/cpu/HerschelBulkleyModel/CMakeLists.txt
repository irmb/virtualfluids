########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(HerschelBulkleyModel)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} FILES hbflow.cpp )