cmake_minimum_required(VERSION 2.6)

PROJECT(CbVectorPoolTests)

#################################################################
# MACHINE_SPECIFIC CMAKE_CONFIG_FILE
#################################################################
INCLUDE("../../../../../../../CMake/CMakeCABMacros.txt")

#################################################################
###   PACKAGES                                               ###
#################################################################
INCLUDE(${SOURCE_ROOT}/basics/container/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/utilities/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/memory/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/objects/CMakePackage.txt)

#################################################################
###   OWN DEFINES 						###
#################################################################
FILE(GLOB SPECIFIC_FILES ${SOURCE_ROOT}/basics/container/examples/CbVectorPool/*.h
                         ${SOURCE_ROOT}/basics/container/examples/CbVectorPool/*.cpp ) 

SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
SOURCE_GROUP(z_specific FILES ${SPECIFIC_FILES})

#################################################################
###  PROJECT ERSTELLEN                                        ###
#################################################################
CREATE_CAB_PROJECT()
