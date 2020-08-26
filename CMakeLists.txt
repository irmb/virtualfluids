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
cmake_minimum_required(VERSION 3.13..3.17 FATAL_ERROR)

project(VirtualFluids CXX)

set (CMAKE_CXX_STANDARD 11)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ".cmake")
set(libraryFolder "libs")
set(testFolder    "tests")
set(appFolder     "apps")
set(thirdFolder   "3rd")

#################################################################################
#  OPTIONS
#################################################################################
option(BUILD_VF_CPU "Build VirtualFluids cpu variant" OFF)
option(BUILD_VF_GPU "Build VirtualFluids gpu variant" OFF)
option(BUILD_VF_UNIT_TESTS "Build VirtualFluids unit tests" OFF)

option(BUILD_SHARED_LIBS "" ON)
#################################################################################
#  MACROS
#################################################################################
include(CMakePrintHelpers)
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
    add_subdirectory(${CMAKE_SOURCE_DIR}/3rdParty/googletest)
endif()
