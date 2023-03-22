PROJECT(RisingBubble2D)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} MultiphaseFlow FILES RisingBubble2D.cpp )
