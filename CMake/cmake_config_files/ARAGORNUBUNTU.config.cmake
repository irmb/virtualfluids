#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Anna Wellmann
# OS:          Ubuntu 20.04 (Docker container)
#################################################################################

set(CMAKE_CUDA_ARCHITECTURES 86)     # Nvidia GeForce RTX 3060

set(GPU_APP "apps/gpu/LBM/")
list(APPEND USER_APPS 
    "${GPU_APP}DrivenCavityMultiGPU"
    "${GPU_APP}SphereScaling"
    # "${GPU_APP}MusselOyster"
    )