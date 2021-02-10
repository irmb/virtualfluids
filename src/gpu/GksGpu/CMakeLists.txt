project(GksGpu LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE shared PRIVATE_LINK basics GksMeshAdapter OpenMP::OpenMP_CXX MPI::MPI_CXX)
