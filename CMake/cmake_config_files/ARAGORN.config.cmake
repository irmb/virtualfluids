#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Anna Wellmann
# OS:          Windows 11
#################################################################################

set(CMAKE_CUDA_ARCHITECTURES 86)     # Nvidia GeForce RTX 3060

# numerical tests location of the grids
# SET(PATH_NUMERICAL_TESTS "E:/temp/numericalTests/")
# list(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")

# add invidual apps here
set(GPU_APP "apps/gpu/LBM/")
list(APPEND USER_APPS 
    "${GPU_APP}DrivenCavityMultiGPU"
    "${GPU_APP}SphereScaling"
    # "${GPU_APP}MusselOyster"
    )
