CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(vfscript)

INCLUDE(${SOURCE_ROOT}/core/IncludsList.txt) 

#################################################################
###   LOCAL FILES                                             ###
#################################################################
FILE(GLOB SPECIFIC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h
                         ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
                         ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp  )
 
SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
SOURCE_GROUP(src FILES ${SPECIFIC_FILES})
  
SET(CAB_ADDITIONAL_LINK_LIBRARIES core)

#YAML support
SET(LINK_LIBRARY optimized ${YAML_RELEASE_LIBRARY} debug ${YAML_DEBUG_LIBRARY})
SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})

#################################################################
###   CREATE PROJECT                                          ###
#################################################################
CREATE_CAB_PROJECT(vfscript BINARY)
