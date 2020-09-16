
if(UNIX)
    set(CMAKE_CXX_STANDARD 14)
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


#############################################################
###                         OPTIONS                       ###
#############################################################
option(VF.BUILD_VF_GPU          "Build VirtualFluids GPU"     ON )
option(VF.BUILD_VF_GKS          "Build VirtualFluids GKS"     OFF )
option(VF.BUILD_VF_TRAFFIC      "Build VirtualFluids Traffic" ON)
option(VF.BUILD_JSONCPP         "Builds json cpp "            OFF)
option(VF.BUILD_NUMERIC_TESTS   "Build numeric tests"         OFF)

option(VF.BUILD_DOUBLE_ACCURACY "Use double accuracy"         OFF )

IF( VF.BUILD_DOUBLE_ACCURACY )
    list(APPEND VF_COMPILER_DEFINITION VF_DOUBLE_ACCURACY)
ENDIF()

#############################################################

enable_language(CUDA)

#############################################################

IF( VF.BUILD_VF_GKS )
    # only use this with device of CC larger than 6.0
    set(CMAKE_CUDA_FLAGS "-Xptxas=\"-v\"" CACHE STRING "" FORCE)
    set(CMAKE_CUDA_ARCHITECTURES 60)
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
    add_subdirectory(src/gpu/GksMeshAdapter)
    add_subdirectory(src/gpu/GksVtkAdapter)

    add_subdirectory(src/gpu/GksGpu)

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
    add_subdirectory(apps/gpu/GKS/Flame7cm)
    #add_subdirectory(targets/apps/GKS/SandiaFlame_1m)
    #add_subdirectory(targets/apps/GKS/Candle)

    #add_subdirectory(targets/apps/GKS/MultiGPU)
    #add_subdirectory(targets/apps/GKS/MultiGPU_nD)
    #add_subdirectory(targets/apps/GKS/SingleGPU)
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