PROJECT(Nozzle)

vf_add_library(BUILDTYPE binary PRIVATE_LINK cpu_core basics ${MPI_CXX_LIBRARIES} LiggghtsCoupling MultiphaseFlow NonNewtonianFluids FILES nozzleSinglePhase.cpp )
