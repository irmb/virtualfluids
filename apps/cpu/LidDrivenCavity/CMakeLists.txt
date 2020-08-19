

PROJECT(LidDrivenCavity)

INCLUDE(${APPS_ROOT_CPU}/IncludsList.cmake)


vf_add_library(BUILDTYPE binary DEPENDS VirtualFluidsCore basics FILES LidDrivenCavity.cpp)
