

SET(USE_INTEL OFF CACHE BOOL "include Intel compiler support")
SET(USE_GCC ON CACHE BOOL "include gcc compiler support")


set (SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(SOURCE_ROOT ${CMAKE_SOURCE_DIR})

#CAB
include("CMake/CMakeCABMacros.cmake") #TODO: Currently we have to include the CABMacros also here, so that the USE_* are defined in the config files for the cpu version

add_subdirectory(3rdParty/MuParser)

add_subdirectory(src/cpu/VirtualFluidsCore)

set (APPS_ROOT_CPU "${CMAKE_SOURCE_DIR}/apps/cpu/")
include(${APPS_ROOT_CPU}/Applications.cmake)