PROJECT(DrivenCavity)

#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-DOMPI_SKIP_MPICXX" )

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluids_GPU basics GridGenerator FILES DrivenCavity.cpp )

include (${CMAKE_PATH}/3rd/cuda.cmake)
include (${CMAKE_PATH}/3rd/mpi.cmake)
include (${CMAKE_PATH}/3rd/boost.cmake)
linkBoost ("serialization")

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/basics)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/gpu)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/gpu/VirtualFluids_GPU)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/gpu/GridGenerator)
target_include_directories(${library_name} PRIVATE ${SOURCE_PATH}/basics/Core)
