
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

option(VF_DOUBLE_ACCURACY       "Use double accuracy"     ON )


#############################################################

enable_language(CUDA)

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
