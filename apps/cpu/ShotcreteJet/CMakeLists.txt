PROJECT(ShotcreteJet)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} LiggghtsCoupling MultiphaseFlow NonNewtonianFluids FILES jet.cpp )
