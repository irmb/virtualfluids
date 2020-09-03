
vf_add_library(BUILDTYPE shared DEPENDS basics GksMeshAdapter)

include (${VF_CMAKE_DIR}/3rd/cuda.cmake)
include (${VF_CMAKE_DIR}/3rd/mpi.cmake)
include (${VF_CMAKE_DIR}/3rd/OpenMP.cmake)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu/GksMeshAdapter")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics/Core")
