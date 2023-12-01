#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Anna Wellmann
# OS:          Windows 11
#################################################################################

set(CMAKE_CUDA_ARCHITECTURES 86)     # Nvidia GeForce RTX 3060

# add invidual apps here
set(GPU_APP "apps/gpu/")
list(APPEND USER_APPS 
    "${GPU_APP}DrivenCavityMultiGPU"
    "${GPU_APP}SphereMultiGPU"
    )
