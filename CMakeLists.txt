cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

if(POLICY CMP0042)
   CMAKE_POLICY(SET CMP0042 NEW)
endif()
if(POLICY CMP0020)
   CMAKE_POLICY(SET CMP0020 NEW)
endif()
if(POLICY CMP0028)
   CMAKE_POLICY(SET CMP0028 NEW)
endif()
if(POLICY CMP0037)
   CMAKE_POLICY(SET CMP0037 NEW)
endif()
if(POLICY CMP0047)
   CMAKE_POLICY(SET CMP0047 NEW)
endif()
if(POLICY CMP0053)
   CMAKE_POLICY(SET CMP0053 NEW)
endif()
if(POLICY CMP0054)
   CMAKE_POLICY(SET CMP0054 NEW)
endif()


if(APPLE)
  set(CMAKE_CXX_STANDARD 11)
endif()

#############################################################
###                     ENVIRONMENT                       ###
#############################################################
set(cmakeMacroPath "CMakeMacros")
include(${cmakeMacroPath}/Environment/environment.cmake)

#############################################################
###                   GENERAL MACROS                      ###
#############################################################
include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/BuildTarget.cmake)
include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/BuildTargetUtilities.cmake)
include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/EndingsToCollect.cmake)
include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/FileUtilities.cmake)

#############################################################
###                   PROJECT SETTINGS                    ###
#############################################################
project(HULC)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_BINARY_DIR})

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ".cmake")
set(libraryFolder "libs")
set(testFolder "tests")
set(appFolder "apps")
set(thirdPartyFolder "3rdParty")

IF(MSVC)
   ADD_DEFINITIONS ( "-DNOMINMAX" )                # Disable Min/Max-Macros
   ADD_DEFINITIONS ( "-D_CRT_SECURE_NO_WARNINGS" ) # disable warnings promoting Microsoft's security enhanced CRT
   ADD_DEFINITIONS ( "-D_SCL_SECURE_NO_WARNINGS" ) # disable warnings triggered by Microsoft's checked iterators
   SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -MP" ) # enable multi-threaded compiling
   SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /bigobj" ) # enable big object files (fatal error C1128)
ENDIF(MSVC)

#############################################################
###                         OPTIONS                       ###
#############################################################
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(HULC.BUILD_VirtualFluids "Builds virtual Fluids" OFF)
option(HULC.BUILD_dem_coupling "Builds dem coupling " ON)
option(HULC.BUILD_GPU_CODE "Build gpu VirtualFluids" ON)
option(HULC.BUILD_JSONCPP "Builds json cpp " ON)
option(HULC.BUILD_TESTS "Build shared libraries" ON)
option(HULC.BUILD_NUMERIC_TESTS "Build shared libraries" ON)

option(VF_DOUBLE_ACCURACY "Build shared libraries" ON)

if(HULC.BUILD_GPU_CODE)
  #SET(CMAKE_CUDA_COMPILER $ENV{CUDA_PATH_V9_0}/bin/nvcc.exe)
  enable_language(CUDA)
endif()
sharedLibs()

include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/FindCompiler.cmake)
configure_file(src/VirtualFluidsDefinitions.in.h VirtualFluidsDefinitions.h)
if(MSVC)
  SET( CMAKE_CXX_FLAGS "/FI${CMAKE_BINARY_DIR}/VirtualFluidsDefinitions.h ${CMAKE_CXX_FLAGS}" )
ELSE(MSVC)
  SET( CMAKE_CXX_FLAGS "-include ${CMAKE_BINARY_DIR}/VirtualFluidsDefinitions.h ${CMAKE_CXX_FLAGS}" )
ENDIF(MSVC)

#############################################################
###                  Virtual Fluids                       ###
#############################################################
IF (HULC.BUILD_VirtualFluids)
  add_definitions(-DVF_METIS)
  add_definitions(-DVF_MPI)
  
  add_subdirectory(3rdParty/MuParser)
  add_subdirectory(targets/libs/VirtualFluids)
ELSE()
  MESSAGE( STATUS "exclude Virtual Fluids." )
ENDIF()


#############################################################
###                dem_coupling                           ###
#############################################################
IF (HULC.BUILD_dem_coupling AND HULC.BUILD_VirtualFluids)
  link_directories("${PE_BINARY_DIR}/src/pe" "${PE_BINARY_DIR}/src/core" "${PE_BINARY_DIR}/src/blockforest" "${PE_BINARY_DIR}/src/domain_decomposition" )
IF(MSVC)
   SET( CMAKE_CXX_FLAGS "/FI ${PE_BINARY_DIR}/src/pe/extern/libccd/ccd/config.h ${CMAKE_CXX_FLAGS}" )
ELSE(MSVC)
   SET( CMAKE_CXX_FLAGS "-include ${PE_BINARY_DIR}/src/pe/extern/libccd/ccd/config.h ${CMAKE_CXX_FLAGS}" )
ENDIF(MSVC)
   
   add_subdirectory(targets/libs/dem_coupling)
ELSE()
  MESSAGE(STATUS "exclude dem-coupling." )
ENDIF()

#############################################################
###                  Virtual Fluids GPU                   ###
#############################################################
IF (HULC.BUILD_GPU_CODE)
 add_subdirectory(targets/libs/VirtualFluids_GPU)
 add_subdirectory(targets/libs/GridGenerator)
ELSE()
  MESSAGE( STATUS "exclude Virtual Fluids GPU." )
ENDIF()

#############################################################
###                     JSONCPP                           ###
#############################################################
IF (NOT HULC.BUILD_JSONCPP)
   MESSAGE( STATUS "Build Input Project without JsonCpp." )
ELSE()
   add_subdirectory(3rdParty/jsoncpp)
   add_definitions(-DBUILD_JSONCPP)
ENDIF()

#############################################################
###                   SUBDIRECTORIES                      ###
#############################################################

add_subdirectory(targets/apps/HULC)
add_subdirectory(targets/libs/VirtualFluidsBasics)
add_subdirectory(targets/libs/core)

if(HULC.BUILD_GPU_CODE OR HULC.BUILD_VirtualFluids)
 add_subdirectory(3rdParty/metis/metis-5.1.0)
endif()


if(HULC.BUILD_TESTS)
  if(VF_DOUBLE_ACCURACY)
    add_definitions("-DRealEq=testing::DoubleEq")
  else()
    add_definitions("-DRealEq=testing::FloatEq")
  endif()

  add_subdirectory(3rdParty/googletest)
  if (HULC.BUILD_GPU_CODE)
    add_subdirectory(targets/tests/GridGeneratorTest)
	add_subdirectory(targets/tests/VirtualFluids_GPUTest)
  endif()
  
  if (HULC.BUILD_VirtualFluids)
    add_subdirectory(targets/tests/VirtualFluidsTest)
  endif()
  
  add_subdirectory(targets/tests/coreTest)
  add_subdirectory(targets/tests/VirtualFluidsBasicsTest)
endif()

if(HULC.BUILD_NUMERIC_TESTS)
 	add_subdirectory(3rdParty/fftw/fftw-3.3.7)
	add_subdirectory(targets/tests/TestingHULC)
endif()

