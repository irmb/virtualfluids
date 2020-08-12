cmake_minimum_required(VERSION 3.9..3.17 FATAL_ERROR)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(VirtualFluids CXX)

set (CMAKE_PATH "${CMAKE_SOURCE_DIR}/CMake")

option(BUILD_VF_CPU "Build VirtualFluids cpu variant" OFF)
option(BUILD_VF_GPU "Build VirtualFluids gpu variant" ON)


include("${CMAKE_PATH}/CMakeCABMacros.cmake")
include("${CMAKE_PATH}/FileUtilities.cmake")
include("${CMAKE_PATH}/VirtualFluidsMacros.cmake")




add_subdirectory(src/basics)

#if (BUILD_VF_CPU)
#include (cpu.cmake)
#endif()
#if(BUILD_VF_GPU)
include (gpu.cmake)
#endif()

