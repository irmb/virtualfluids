
vf_add_library(BUILDTYPE static DEPENDS basics GksMeshAdapter)

include (${VF_CMAKE_DIR}/3rd/cuda.cmake)
include (${VF_CMAKE_DIR}/3rd/mpi.cmake)
include (${VF_CMAKE_DIR}/3rd/OpenMP.cmake)
