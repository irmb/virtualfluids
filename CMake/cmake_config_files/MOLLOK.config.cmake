#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Martin Schoenherr
# OS:          Windows 10
#################################################################################

#SET TO CORRECT PATH:
SET(CMAKE_CUDA_ARCHITECTURES 86)

SET(PATH_NUMERICAL_TESTS "D:/out/numericalTests/")
LIST(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")
