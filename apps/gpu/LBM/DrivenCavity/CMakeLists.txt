PROJECT(DrivenCavity)

#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-DOMPI_SKIP_MPICXX" )

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluids_GPU basics GridGenerator FILES DrivenCavity.cpp )

include (${VF_CMAKE_DIR}/3rd/cuda.cmake)
include (${VF_CMAKE_DIR}/3rd/mpi.cmake)
include (${VF_CMAKE_DIR}/3rd/boost.cmake)
linkBoost ("serialization")

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/gpu)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/gpu/VirtualFluids_GPU)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/gpu/GridGenerator)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/Core)
