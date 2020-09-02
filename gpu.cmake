cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

if(POLICY CMP0042)
    CMAKE_POLICY(SET CMP0042 NEW)
endif()
if(POLICY CMP0020)
    CMAKE_POLICY(SET CMP0020 NEW)
endif()
if(POLICY CMP0028)
    CMAKE_POLICY(SET CMP0028 NEW)
endif()
if(POLICY CMP0037)
    CMAKE_POLICY(SET CMP0037 NEW)
endif()
if(POLICY CMP0047)
    CMAKE_POLICY(SET CMP0047 NEW)
endif()
if(POLICY CMP0053)
    CMAKE_POLICY(SET CMP0053 NEW)
endif()
if(POLICY CMP0054)
    CMAKE_POLICY(SET CMP0054 NEW)
endif()


if(UNIX)
    set(CMAKE_CXX_STANDARD 11)
endif()

#############################################################
###                     CUDAPATH                          ###
#############################################################

# if CMake cannot find CUDA by itself, set the correct paths manually:
#SET(CUDA_CUT_INCLUDE_DIR    "/cluster/cuda/9.0/include;/cluster/cuda/9.0/samples/common/inc" CACHE PATH "CUDA_CUT_INCLUDE_DIR")
#SET(CUDA_SAMPLE_INCLUDE_DIR "/cluster/cuda/9.0/samples/common/inc" CACHE PATH "CUDA_CUT_INCLUDE_DIR")

#############################################################
###                   PROJECT SETTINGS                    ###
#############################################################

project(VirtualFluidsGPU)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_BINARY_DIR}/gpu)
#
#set(libraryFolder    "libs")
#set(gksLibraryFolder "libs/GKS")
#
#set(testFolder "tests")
#
#set(appFolder    "apps")
#set(lbmAppFolder "apps/LBM")
#set(gksAppFolder "apps/GKS")
#
#set(thirdPartyFolder "3rdParty")

IF(MSVC)
    ADD_DEFINITIONS ( "-DNOMINMAX" )                # Disable Min/Max-Macros
    ADD_DEFINITIONS ( "-D_CRT_SECURE_NO_WARNINGS" ) # disable warnings promoting Microsoft's security enhanced CRT
    ADD_DEFINITIONS ( "-D_SCL_SECURE_NO_WARNINGS" ) # disable warnings triggered by Microsoft's checked iterators
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -MP" ) # enable multi-threaded compiling
    SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /bigobj" ) # enable big object files (fatal error C1128)
ENDIF(MSVC)

#############################################################
###                         OPTIONS                       ###
#############################################################
option(BUILD_SHARED_LIBS        "Build shared libraries"  ON )

option(VF_DOUBLE_ACCURACY       "Use double accuracy"     ON )


#############################################################

enable_language(CUDA)

#sharedLibs()

#############################################################


# only use this with device of CC larger than 6.0
IF(VF_DOUBLE_ACCURACY)
    set(CMAKE_CUDA_FLAGS " -arch=sm_60" CACHE STRING "" FORCE)
ELSE(VF_DOUBLE_ACCURACY)
    set(CMAKE_CUDA_FLAGS "" CACHE STRING "" FORCE)
ENDIF(VF_DOUBLE_ACCURACY)
set(CMAKE_CUDA_FLAGS_DEBUG " -G" CACHE STRING "" FORCE)

#############################################################
###                  Virtual Fluids GPU                   ###
#############################################################

add_subdirectory(src/gpu/GridGenerator)
add_subdirectory(src/gpu/VirtualFluids_GPU)

add_subdirectory(src/gpu/GksMeshAdapter)
add_subdirectory(src/gpu/GksGpu)

add_subdirectory(apps/gpu/LidDrivenCavity)
