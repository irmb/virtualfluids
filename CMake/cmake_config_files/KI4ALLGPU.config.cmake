#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Anna Wellmann
# OS:          Ubuntu 20.04 (Docker container)
#################################################################################

set(CMAKE_CUDA_ARCHITECTURES 80)     # Nvidia Tesla A100

set(GPU_APP "apps/gpu/")
list(APPEND USER_APPS 
    "${GPU_APP}SphereMultiGPU"
)

# run docker container with:
# docker run -it -v `pwd`:`pwd` -w `pwd` --gpus all --hostname ki4allgpu --name virtual-fluids-environment <containerid>