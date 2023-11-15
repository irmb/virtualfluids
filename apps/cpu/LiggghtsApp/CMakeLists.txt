PROJECT(LiggghtsApp)

vf_add_library(BUILDTYPE binary PRIVATE_LINK cpu_core basics ${MPI_CXX_LIBRARIES} LiggghtsCoupling FILES LiggghtsApp.cpp )
