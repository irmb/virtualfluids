PROJECT(FallingSphere)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} LiggghtsCoupling FILES FallingSphere.cpp )
