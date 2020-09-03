PROJECT(PChannel)

INCLUDE(${APPS_ROOT}/IncludsList.cmake) 

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluidsCore VirtualFluidsBasic FILES pChannel.cpp)
