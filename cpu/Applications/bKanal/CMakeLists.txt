CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(bKanal)

INCLUDE(${SOURCE_ROOT}/lib/IncludsList.txt) 

#################################################################
###   LOCAL FILES                                             ###
#################################################################
FILE(GLOB SPECIFIC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h
                         ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
                         ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp  )
 
SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
SOURCE_GROUP(src FILES ${SPECIFIC_FILES})
  
SET(CAB_ADDITIONAL_LINK_LIBRARIES vfluids)

#################################################################
###   CREATE PROJECT                                          ###
#################################################################
CREATE_CAB_PROJECT(bKanal BINARY)
