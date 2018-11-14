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
option(BUILD_SHARED_LIBS        "Build shared libraries"  ON )
option(HULC.BUILD_GPU_CODE      "Build gpu VirtualFluids" ON )
option(HULC.BUILD_JSONCPP       "Builds json cpp "        ON )
option(HULC.BUILD_NUMERIC_TESTS "Build numeric tests"     OFF)

option(VF_DOUBLE_ACCURACY       "Build shared libraries"  OFF)

#############################################################

if(HULC.BUILD_GPU_CODE)
  #SET(CMAKE_CUDA_COMPILER $ENV{CUDA_PATH_V9_0}/bin/nvcc.exe)
  enable_language(CUDA)
endif()
sharedLibs()

#############################################################

include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/FindCompiler.cmake)
configure_file(src/VirtualFluidsDefinitions.in.h VirtualFluidsDefinitions.h)
if(MSVC)
  SET( CMAKE_CXX_FLAGS "/FI${CMAKE_BINARY_DIR}/VirtualFluidsDefinitions.h ${CMAKE_CXX_FLAGS}" )
ELSE(MSVC)
  SET( CMAKE_CXX_FLAGS "-include ${CMAKE_BINARY_DIR}/VirtualFluidsDefinitions.h ${CMAKE_CXX_FLAGS}" )
ENDIF(MSVC)

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

add_subdirectory(3rdParty/metis/metis-5.1.0)

add_subdirectory(targets/apps/HULC)
add_subdirectory(targets/libs/VirtualFluidsBasics)
add_subdirectory(targets/libs/Core)

if(HULC.BUILD_NUMERIC_TESTS)
 	add_subdirectory(3rdParty/fftw/fftw-3.3.7)
	add_subdirectory(targets/tests/TestingHULC)
endif()

