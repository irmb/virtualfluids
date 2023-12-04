#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Anna Wellmann
# OS:          Ubuntu 22.04
#################################################################################

set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
SET(CMAKE_CUDA_ARCHITECTURES 61)    # GeForce 1080 Ti

set(GPU_APP "apps/gpu/")
list(APPEND USER_APPS 
    "${GPU_APP}DrivenCavityMultiGPU"
    "${GPU_APP}SphereMultiGPU"
    )