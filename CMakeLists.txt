#################################################################################
#   _    ___      __              __________      _     __
# | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
# | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
# | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
# |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
#
#################################################################################
#  required cmake versions
#  CMAKE 3.9: CUDA support
#################################################################################
cmake_minimum_required(VERSION 3.9..3.17 FATAL_ERROR)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(VirtualFluids CXX)

#################################################################################
#  OPTIONS
#################################################################################
option(BUILD_VF_CPU "Build VirtualFluids cpu variant" OFF)
option(BUILD_VF_GPU "Build VirtualFluids gpu variant" ON)
option(BUILD_VF_UNIT_TESTS "Build VirtualFluids unit tests" ON)

#################################################################################
#  MACROS
#################################################################################
set (CMAKE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include("${CMAKE_PATH}/CMakeCABMacros.cmake")
include("${CMAKE_PATH}/FileUtilities.cmake")
include("${CMAKE_PATH}/VirtualFluidsMacros.cmake")

#################################################################################
#  COMMON LIBRARIES
#################################################################################
add_subdirectory(src/basics)

#################################################################################
#  VIRTUAL FLUIDS CPU / GPU
#################################################################################
#if (BUILD_VF_CPU)
include (cpu.cmake)
#endif()
#if(BUILD_VF_GPU)
#include (gpu.cmake)
#endif()

