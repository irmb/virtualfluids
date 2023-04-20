########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(HerschelBulkleySphere)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} NonNewtonianFluids FILES hbsphere.cpp )