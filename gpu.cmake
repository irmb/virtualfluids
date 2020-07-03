cmake_minimum_required(VERSION 3.9..3.17 FATAL_ERROR)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(VirtualFluids CXX)

set (CMAKE_PATH "${CMAKE_SOURCE_DIR}/CMake")

option(BUILD_VF_CPU "Build VirtualFluids cpu variant" OFF)
option(BUILD_VF_GPU "Build VirtualFluids gpu variant" ON)


include("${CMAKE_PATH}/CMakeCABMacros.cmake")
include("${CMAKE_PATH}/FileUtilities.cmake")
include("${CMAKE_PATH}/VirtualFluidsMacros.cmake")


SET(USE_ZOLTAN OFF CACHE BOOL "include Zoltan library support")
SET(USE_METIS ON CACHE BOOL "include METIS library support")
SET(USE_MPI ON CACHE BOOL "include MPI library support")
SET(USE_VTK OFF CACHE BOOL "include VTK library support")
SET(USE_CATALYST OFF CACHE BOOL "include Paraview Catalyst support")
SET(USE_BOOST OFF CACHE BOOL "include Boost support")
#SET(USE_PYTHON OFF CACHE BOOL "include Python scripting support")
#SET(USE_FETOL OFF CACHE BOOL "include FETOL library support")
SET(USE_INTEL OFF CACHE BOOL "include Intel compiler support")
SET(USE_GCC OFF CACHE BOOL "include gcc compiler support")
SET(USE_HLRN_LUSTRE OFF CACHE BOOL "include HLRN Lustre support")
SET(USE_DEM_COUPLING OFF CACHE BOOL "PE plugin")
IF(${USE_MPI})
    LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_MPI)
ENDIF()
IF(${USE_METIS})
    LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_METIS)
ENDIF()
# FIND MPI
IF((NOT ${CMAKE_CXX_COMPILER} MATCHES mpicxx) AND (NOT ${CMAKE_CXX_COMPILER} MATCHES mpiicpc))# OR NOT ${CMAKE_CXX_COMPILER} MATCHES cc OR NOT ${CMAKE_CXX_COMPILER} MATCHES mpiCC)
    FIND_PACKAGE(MPI REQUIRED)
ENDIF()

add_subdirectory(src/basics)

#if (BUILD_VF_CPU)
    add_subdirectory(3rdParty/MuParser)
    add_subdirectory(cpu)
#endif()
#if(BUILD_VF_GPU)
#    add_subdirectory(gpu)
#endif()

set (APPS_ROOT_CPU "${CMAKE_SOURCE_DIR}/apps/cpu/")
include(${APPS_ROOT_CPU}/Applications.cmake)