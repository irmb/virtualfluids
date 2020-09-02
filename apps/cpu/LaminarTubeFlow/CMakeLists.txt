
PROJECT(ltf)

INCLUDE(${APPS_ROOT_CPU}/IncludsList.cmake)


vf_add_library(BUILDTYPE binary DEPENDS VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} FILES ltf.cpp )
