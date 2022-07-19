#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Martin Schoenherr
# OS:          Windows 10
#################################################################################

#SET TO CORRECT PATH:
set(CMAKE_CUDA_ARCHITECTURES 86)

set(PATH_NUMERICAL_TESTS "D:/out/numericalTests/")
list(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")


list(APPEND USER_APPS "apps/gpu/LBM/WTG_RUB")