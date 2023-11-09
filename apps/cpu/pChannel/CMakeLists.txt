PROJECT(PChannel)

INCLUDE(${APPS_ROOT}/IncludsList.cmake) 

vf_add_library(BUILDTYPE binary PRIVATE_LINK cpu_core VirtualFluidsBasic FILES pChannel.cpp)
