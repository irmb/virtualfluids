PROJECT(Nozzle)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} LiggghtsCoupling FILES nozzle.cpp )
