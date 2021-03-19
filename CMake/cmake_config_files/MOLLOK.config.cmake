#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Martin Schoenherr
# OS:          Windows 10
#################################################################################

#SET TO CORRECT PATH:
SET(BOOST_ROOT  "D:/libraries/boost_1_74_0"  CACHE PATH "BOOST_ROOT")
SET(BOOST_LIBRARYDIR  "D:/libraries/boost_1_74_0/stageMSVC64VS2019/lib" CACHE PATH "BOOST_LIBRARYDIR")
SET(CMAKE_CUDA_ARCHITECTURES 52)

SET(PATH_NUMERICAL_TESTS "E:/temp/numericalTests/")
LIST(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")
