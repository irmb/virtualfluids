PROJECT(GridGeneratorTest)

vf_add_library(BUILDTYPE binary DEPENDS basics VirtualFluids_GPU GridGenerator FILES gridGenerator.cpp)

include (${VF_CMAKE_DIR}/3rd/cuda.cmake)
include (${VF_CMAKE_DIR}/3rd/mpi.cmake)
include (${VF_CMAKE_DIR}/3rd/boost.cmake)
linkBoost ("serialization")

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu/VirtualFluids_GPU")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu/GridGenerator")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics/Core")