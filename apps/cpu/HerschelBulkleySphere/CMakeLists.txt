CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(HerschelBulkleySphere)

INCLUDE(${APPS_ROOT_CPU}/IncludsList.cmake) 

################################################################
##   LOCAL FILES                                             ###
################################################################
# FILE(GLOB SPECIFIC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h
                         # ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
                         # ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp  )
 
# SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
# SOURCE_GROUP(src FILES ${SPECIFIC_FILES})
  
# SET(CAB_ADDITIONAL_LINK_LIBRARIES VirtualFluids)

################################################################
##   CREATE PROJECT                                          ###
################################################################
# CREATE_CAB_PROJECT(hbsphere BINARY)

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluidsCore basics ${MPI_CXX_LIBRARIES} FILES hbsphere.cpp )