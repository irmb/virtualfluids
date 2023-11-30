#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Konstantin Kutscher
# OS:          CentOS 7.3
#################################################################################


## nvidia
set(CMAKE_CUDA_ARCHITECTURES 60) # NVIDIA Tesla P100

set(GPU_APP "apps/gpu/")
list(APPEND USER_APPS 
    # "${GPU_APP}DrivenCavityMultiGPU"
    # "${GPU_APP}SphereMultiGPU"
    # "${GPU_APP}MusselOyster"
    )
