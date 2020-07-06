

#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Cuda/Link.cmake)
#linkCuda(${targetName})
#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/MPI/Link.cmake)
#linkMPI(${targetName})
#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Boost/Link.cmake)
#linkBoost(${targetName} "serialization")
#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Metis/Link.cmake)
#linkMetis(${targetName})
#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/OpenMP/Link.cmake)
#linkOpenMP(${targetName})


vf_add_library(BUILDTYPE shared)

include (${CMAKE_PATH}/3rd/cuda.cmake)
include (${CMAKE_PATH}/3rd/OpenMP.cmake)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics/Core")


set_property(TARGET ${library_name} PROPERTY CUDA_SEPARABLE_COMPILATION ON)