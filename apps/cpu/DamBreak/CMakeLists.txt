PROJECT(DamBreak)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} MultiphaseFlow NonNewtonianFluids FILES main.cpp )
