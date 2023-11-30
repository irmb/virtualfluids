#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Martin Schoenherr
# OS:          Windows 11
#################################################################################

# cuda compute capability
set(CMAKE_CUDA_ARCHITECTURES 86)

# numerical tests location of the grids
set(PATH_NUMERICAL_TESTS "D:/out/numericalTests/")
list(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")
