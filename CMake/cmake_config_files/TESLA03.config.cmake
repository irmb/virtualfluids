#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Martin Schoenherr
# OS:          Windows 10
#################################################################################

# cuda compute capability
SET(CMAKE_CUDA_ARCHITECTURES 52)

# numerical tests location of the grids
SET(PATH_NUMERICAL_TESTS "E:/temp/numericalTests/")
list(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")
