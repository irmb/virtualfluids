#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Anna Wellmann
# OS:          Ubuntu 20.04 (Docker container)
#################################################################################

set(CMAKE_CUDA_ARCHITECTURES 86)     # Nvidia GeForce RTX 3060

set(PATH_NUMERICAL_TESTS "D:/out/numericalTests/")
list(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")

set(GPU_APP "apps/gpu/")
list(APPEND USER_APPS 
    "${GPU_APP}DrivenCavityMultiGPU"
    "${GPU_APP}SphereScaling"
    # "${GPU_APP}MusselOyster"
    )

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)
# add_compile_options(-fsanitize=undefined)
# add_link_options(-fsanitize=undefined)
