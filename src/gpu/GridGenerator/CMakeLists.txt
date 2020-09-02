

vf_add_library(BUILDTYPE shared DEPENDS basics)

include (${CMAKE_PATH}/3rd/cuda.cmake)
include (${CMAKE_PATH}/3rd/OpenMP.cmake)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/basics)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/basics/geometry3d)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/gpu)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/basics/Core)


set_property(TARGET ${library_name} PROPERTY CUDA_SEPARABLE_COMPILATION ON)