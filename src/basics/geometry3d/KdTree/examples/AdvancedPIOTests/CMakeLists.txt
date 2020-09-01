cmake_minimum_required(VERSION 2.6)

PROJECT(kdtree_pioTest)

#################################################################
# MACHINE_SPECIFIC CMAKE_CONFIG_FILE
#################################################################
INCLUDE("../../../../../CMake/CMakeCABMacros.txt")

#################################################################
###   PACKAGES                                               ###
#################################################################
INCLUDE(${SOURCE_ROOT}/basics/utilities/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/3rdParty/dirstream/CMakePackage.txt)

INCLUDE(${SOURCE_ROOT}/basics/utilities/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/container/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/memory/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/objects/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/writer/CMakePackage.txt)

INCLUDE(${SOURCE_ROOT}/geometry3d/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/geometry3d/KdTree/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/geometry3d/fem/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/geometry3d/creator/CMakePackage.txt)

#################################################################
###   OWN DEFINES 						###
#################################################################
FILE(GLOB SPECIFIC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h
                         ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ) 

SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
SOURCE_GROUP(z_main FILES ${SPECIFIC_FILES})

#################################################################
###  PROJECT ERSTELLEN                                        ###
#################################################################
CREATE_CAB_PROJECT(${PROJECT_NAME} BINARY)
