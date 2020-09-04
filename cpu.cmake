

SET(USE_INTEL OFF CACHE BOOL "include Intel compiler support")
SET(USE_GCC OFF CACHE BOOL "include gcc compiler support")


#CAB
include("CMake/CMakeCABMacros.cmake") #TODO: Currently we have to include the CABMacros also here, so that the USE_* are defined in the config files for the cpu version

add_subdirectory(${VF_THIRD_DIR}/MuParser)

add_subdirectory(src/cpu/VirtualFluidsCore)

set (APPS_ROOT_CPU "${CMAKE_SOURCE_DIR}/apps/cpu/")
include(${APPS_ROOT_CPU}/Applications.cmake)