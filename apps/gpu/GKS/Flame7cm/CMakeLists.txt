PROJECT(Flame7cm LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics GridGenerator GksMeshAdapter GksVtkAdapter GksGpu MPI::MPI_CXX FILES Flame7cm.cpp )

set_source_files_properties(Flame7cm.cpp PROPERTIES LANGUAGE CUDA)
