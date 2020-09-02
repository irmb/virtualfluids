
PROJECT(LidDrivenCavity)

INCLUDE(${APPS_ROOT}/IncludsList.cmake)

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluidsCore VirtualFluidsBasic FILES LidDrivenCavity.cpp)
