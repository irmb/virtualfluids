#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Martin Schoenherr
# OS:          Windows 10
#################################################################################

#Don't change:
set(METIS_ROOT ${CMAKE_SOURCE_DIR}/3rdParty/metis/metis-5.1.0 CACHE PATH "METIS ROOT")
set(GMOCK_ROOT ${CMAKE_SOURCE_DIR}/3rdParty/googletest CACHE PATH "GMOCK ROOT")
set(JSONCPP_ROOT ${CMAKE_SOURCE_DIR}/3rdParty/jsoncpp CACHE PATH "JSONCPP ROOT")
set(FFTW_ROOT ${CMAKE_SOURCE_DIR}/3rdParty/fftw/fftw-3.3.7 CACHE PATH "JSONCPP ROOT")

set(CMAKE_CUDA_ARCHITECTURES 52)

set(VTK_DIR "F:/Libraries/vtk/VTK-8.2.0/build" CACHE PATH "VTK directory override" FORCE)

set(PATH_NUMERICAL_TESTS "E:/temp/numericalTests/")
list(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")
set(CMAKE_CUDA_ARCHITECTURES 52)
set(CMAKE_CUDA_ARCHITECTURES 52)

list(APPEND USER_APPS "apps/gpu/LBM/MusselOyster")
