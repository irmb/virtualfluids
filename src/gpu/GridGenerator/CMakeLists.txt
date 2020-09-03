

vf_add_library(BUILDTYPE shared DEPENDS basics)

include (${VF_CMAKE_DIR}/3rd/cuda.cmake)
include (${VF_CMAKE_DIR}/3rd/OpenMP.cmake)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/geometry3d)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/gpu)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/Core)


set_property(TARGET ${library_name} PROPERTY CUDA_SEPARABLE_COMPILATION ON)