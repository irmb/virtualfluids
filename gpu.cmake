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
    set(CMAKE_CXX_STANDARD 14)
endif()

#############################################################
###                     ENVIRONMENT                       ###
#############################################################
set(cmakeMacroPath "CMakeMacros")
#include(${cmakeMacroPath}/Environment/environment.cmake)

#CAB
#include("../CMake/CMakeCABMacros.cmake")
#include("../CMake/FileUtilities.cmake")
#include("../CMake/VirtualFluidsMacros.cmake")

#############################################################
###                   GENERAL MACROS                      ###
#############################################################
#include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/BuildTarget.cmake)
#include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/BuildTargetUtilities.cmake)
#include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/EndingsToCollect.cmake)
#include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/FileUtilities.cmake)

#############################################################
###                   PROJECT SETTINGS                    ###
#############################################################

project(VirtualFluidsGPU)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_BINARY_DIR}/gpu)

set(libraryFolder    "libs")
set(gksLibraryFolder "libs/GKS")

set(testFolder "tests")

set(appFolder    "apps")
set(lbmAppFolder "apps/LBM")
set(gksAppFolder "apps/GKS")

set(thirdPartyFolder "3rdParty")

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
option(BUILD_SHARED_LIBS        "Build shared libraries"      ON )
option(VF.BUILD_VF_GPU          "Build VirtualFluids GPU"     ON )
option(VF.BUILD_VF_GKS          "Build VirtualFluids GKS"     OFF )
option(VF.BUILD_VF_TRAFFIC      "Build VirtualFluids Traffic" ON)
option(VF.BUILD_JSONCPP         "Builds json cpp "            OFF)
option(VF.BUILD_NUMERIC_TESTS   "Build numeric tests"         OFF)

option(VF.BUILD_DOUBLE_ACCURACY "Use double accuracy"         OFF )

IF( VF.BUILD_DOUBLE_ACCURACY )
    SET( VF_DOUBLE_ACCURACY 1 )
ENDIF()

#############################################################

enable_language(CUDA)

#sharedLibs()

#############################################################

include(${CMAKE_PATH}/CMakeMacros_old/general/FindCompiler.cmake)
configure_file(src/gpu/VirtualFluidsDefinitions.in.h VirtualFluidsDefinitions.h)
if(MSVC)
    SET( CMAKE_CXX_FLAGS "/FI${CMAKE_BINARY_DIR}/VirtualFluidsDefinitions.h ${CMAKE_CXX_FLAGS}" )
ELSE(MSVC)
    SET( CMAKE_CXX_FLAGS "-include ${CMAKE_BINARY_DIR}/VirtualFluidsDefinitions.h ${CMAKE_CXX_FLAGS}" )
ENDIF(MSVC)

IF( VF.BUILD_VF_GKS )
    # only use this with device of CC larger than 6.0
    set(CMAKE_CUDA_FLAGS " -arch=sm_60 -Xptxas=\"-v\"" CACHE STRING "" FORCE)
ENDIF()

set(CMAKE_CUDA_FLAGS_DEBUG " -G" CACHE STRING "" FORCE)


##########################################################################################################################
###                  Subdirectories                                                                                    ###
##########################################################################################################################

#############################################################
###                  Core                                 ###
#############################################################

add_subdirectory(src/gpu/GridGenerator)
#add_subdirectory(3rdParty/metis/metis-5.1.0)

#############################################################
###                  Virtual Fluids GPU                   ###
#############################################################

IF (VF.BUILD_VF_GPU)
    add_subdirectory(src/gpu/VirtualFluids_GPU)

    #add_subdirectory(targets/apps/LBM/lbmTest)
    #add_subdirectory(targets/apps/LBM/metisTest)
    #add_subdirectory(targets/apps/LBM/Basel)
    #add_subdirectory(targets/apps/LBM/BaselNU)
    #add_subdirectory(targets/apps/LBM/BaselMultiGPU)

    add_subdirectory(apps/gpu/LBM/DrivenCavity)
    add_subdirectory(apps/gpu/LBM/gridGeneratorTest)
    add_subdirectory(apps/gpu/LBM/TGV_3D)
    add_subdirectory(apps/gpu/LBM/TGV_3D_MultiGPU)
ELSE()
    MESSAGE( STATUS "exclude Virtual Fluids GPU." )
ENDIF()

#############################################################
###                  Virtual Fluids GKS                   ###
#############################################################

IF (VF.BUILD_VF_GKS)
    add_subdirectory(targets/libs/GksMeshAdapter)
    add_subdirectory(targets/libs/GksVtkAdapter)

    add_subdirectory(targets/libs/GksGpu)

    #add_subdirectory(targets/apps/GKS/gksTest)
    #add_subdirectory(targets/apps/GKS/ChannelFlow)

    #add_subdirectory(targets/apps/GKS/ChannelFlowObstacle)
    #add_subdirectory(targets/apps/GKS/ShearWave)

    #add_subdirectory(targets/apps/GKS/LiFuXu)

    #add_subdirectory(targets/apps/GKS/TaylorGreen3D)
    #add_subdirectory(targets/apps/GKS/DrivenCavity3D)
    #add_subdirectory(targets/apps/GKS/ThermalCavity)

    #add_subdirectory(targets/apps/GKS/ThermalCavityMultiGPU)
    #add_subdirectory(targets/apps/GKS/DrivenCavityMultiGPU)
    #add_subdirectory(targets/apps/GKS/RayleighBenardMultiGPU)

    #add_subdirectory(targets/apps/GKS/SalinasVazquez)
    #add_subdirectory(targets/apps/GKS/BoundaryJet)

    #add_subdirectory(targets/apps/GKS/PropaneFlame)
    #add_subdirectory(targets/apps/GKS/ConfinedCombustion)
    #add_subdirectory(targets/apps/GKS/MethaneFlame)

    #add_subdirectory(targets/apps/GKS/Room)
    #add_subdirectory(targets/apps/GKS/RoomMultiGPU)
    #add_subdirectory(targets/apps/GKS/RoomFire)
    #add_subdirectory(targets/apps/GKS/RoomFireExtended)
    #add_subdirectory(targets/apps/GKS/ConcreteHeatFluxBCTest)

    #add_subdirectory(targets/apps/GKS/PoolFire)
    add_subdirectory(targets/apps/GKS/Flame7cm)
    add_subdirectory(targets/apps/GKS/SandiaFlame_1m)
    #add_subdirectory(targets/apps/GKS/Candle)

    #add_subdirectory(targets/apps/GKS/MultiGPU)
    add_subdirectory(targets/apps/GKS/MultiGPU_nD)
    add_subdirectory(targets/apps/GKS/SingleGPU)
ELSE()
    MESSAGE( STATUS "exclude Virtual Fluids GKS." )
ENDIF()

#############################################################
###                     JSONCPP                           ###
#############################################################
IF (NOT VF.BUILD_JSONCPP)
    MESSAGE( STATUS "Build Input Project without JsonCpp." )
ELSE()
    add_subdirectory(3rdParty/jsoncpp)
    add_definitions(-DBUILD_JSONCPP)
ENDIF()

#############################################################
###                   Numeric Tests                       ###
#############################################################

if(VF.BUILD_NUMERIC_TESTS)
    add_subdirectory(3rdParty/fftw/fftw-3.3.7)
    add_subdirectory(3rdParty/googletest)
    add_subdirectory(targets/tests/NumericalTests)
    add_subdirectory(targets/tests/NumericalTestPostProcessing)
endif()

#############################################################
###					Annas Traffic Sim				      ###
#############################################################

if(VF.BUILD_VF_TRAFFIC)
    add_subdirectory(src/gpu/Traffic)
    add_subdirectory(apps/gpu/LBM/TrafficTest)
endif()