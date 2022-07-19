#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Soeren Peters
# OS:          Ubuntu 20.04
#################################################################################


SET(CMAKE_CUDA_ARCHITECTURES 70)

list(APPEND USER_APPS "apps/gpu/LBM/WTG_RUB")