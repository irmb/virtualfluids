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
cmake_minimum_required(VERSION 3.13..3.18 FATAL_ERROR)

project(VirtualFluids CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ".cmake")
set(libraryFolder "libs")
set(testFolder    "tests")
set(appFolder     "apps")
set(thirdFolder   "3rd")

set (VF_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set (VF_THIRD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty)
set (VF_SRC_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/src)
set (VF_ROOT_DIR  ${CMAKE_CURRENT_SOURCE_DIR})

#################################################################################
#  OPTIONS
#################################################################################
option(BUILD_VF_CPU "Build VirtualFluids cpu variant" OFF)
option(BUILD_VF_GPU "Build VirtualFluids gpu variant" OFF)
option(BUILD_VF_UNIT_TESTS "Build VirtualFluids unit tests" OFF)

option(BUILD_SHARED_LIBS "" ON)

#################################################################################
#  CMAKE POLICIES
#################################################################################
# CMAKE_CUDA_ARCHITECTURES
# https://cmake.org/cmake/help/git-stage/policy/CMP0104.htmls
if(POLICY CMP0104)
    cmake_policy(SET CMP0104 NEW)
    set(CMAKE_CUDA_ARCHITECTURES 30)
endif()

#################################################################################
#  MACROS
#################################################################################
include(CMakePrintHelpers)
include(${VF_CMAKE_DIR}/VirtualFluidsMacros.cmake)

#################################################################################
#  COMMON LIBRARIES
#################################################################################
add_subdirectory(src/basics)

#################################################################################
#  VIRTUAL FLUIDS CPU / GPU
#################################################################################
if (BUILD_VF_CPU)
    include (cpu.cmake)
endif()
if(BUILD_VF_GPU)
    include (gpu.cmake)
endif()

#################################################################################
#  3rd Party Libraries
#################################################################################
if(BUILD_VF_UNIT_TESTS)
    add_subdirectory(${VF_THIRD_DIR}/googletest)
endif()
