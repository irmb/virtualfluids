CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(cylinder)

INCLUDE(${APPS_ROOT_CPU}/IncludsList.cmake) 

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} FILES cylinder.cpp )
