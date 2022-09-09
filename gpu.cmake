IF( BUILD_VF_GKS )
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

IF (BUILD_VF_GPU)
    add_subdirectory(src/gpu/VirtualFluids_GPU)

    add_subdirectory(apps/gpu/LBM/DrivenCavity)
    add_subdirectory(apps/gpu/LBM/SphereGPU)
    add_subdirectory(apps/gpu/LBM/BoundaryLayer)
ELSE()
    MESSAGE( STATUS "exclude Virtual Fluids GPU." )
ENDIF()

#############################################################
###                  Virtual Fluids GKS                   ###
#############################################################

IF (BUILD_VF_GKS)
    add_subdirectory(src/gpu/GksMeshAdapter)
    add_subdirectory(src/gpu/GksVtkAdapter)

    add_subdirectory(src/gpu/GksGpu)

    add_subdirectory(apps/gpu/GKS/Flame7cm)
ELSE()
    MESSAGE( STATUS "exclude Virtual Fluids GKS." )
ENDIF()

#############################################################
###                     JSONCPP                           ###
#############################################################
IF (NOT BUILD_JSONCPP)
    MESSAGE( STATUS "Build Input Project without JsonCpp." )
ELSE()
    add_subdirectory(3rdParty/jsoncpp)
    add_definitions(-DBUILD_JSONCPP)
ENDIF()

#############################################################
###                   Numeric Tests                       ###
#############################################################

if(BUILD_NUMERIC_TESTS)

    # PATH_NUMERICAL_TESTS can be passed to cmake e.g. cmake .. -DPATH_NUMERICAL_TESTS=/data/
    if(PATH_NUMERICAL_TESTS)
        LIST(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")
    endif()

    if(NOT BUILD_VF_UNIT_TESTS) # in this case googletest is already included.
        add_subdirectory(${VF_THIRD_DIR}/googletest)
    endif()

    add_subdirectory(3rdParty/fftw/fftw-3.3.7)
    add_subdirectory(apps/gpu/tests/NumericalTests)
    add_subdirectory(apps/gpu/tests/NumericalTestPostProcessing)
endif()

#############################################################
###                 Annas Traffic Sim                     ###
#############################################################

if(BUILD_VF_TRAFFIC)
    add_subdirectory(src/gpu/Traffic)
    add_subdirectory(apps/gpu/LBM/TrafficTest)
endif()
