PROJECT(Flame7cm)

vf_add_library(BUILDTYPE binary DEPENDS basics GridGenerator GksMeshAdapter GksVtkAdapter GksGpu FILES Flame7cm.cpp )

include (${VF_CMAKE_DIR}/3rd/cuda.cmake)
include (${VF_CMAKE_DIR}/3rd/mpi.cmake)
