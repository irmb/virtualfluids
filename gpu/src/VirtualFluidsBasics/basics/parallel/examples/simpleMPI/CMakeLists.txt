cmake_minimum_required(VERSION 2.6)

INCLUDE("../../../../CMakeCABMacros.txt")
INCLUDE("../../../../CMakeSetCompilerFlags.txt")
  
CHECK_FOR_VARIABLE(CAB_MACHINE "machine name, e.g. ALTIX, ARWEN")
SET(CMAKE_CONFIG_FILE "${SOURCE_ROOT}/cmake_config_files/${CAB_MACHINE}.config.cmake")

PROJECT(simpleMPI)

#erst hier das config file einfügen, ansonsten werden manche settings durch (Project) überschrieben)  
INCLUDE(${CMAKE_CONFIG_FILE})  
  
SET(EXECUTABLE_NAME simpleMPI)

################################################################
##   PACKAGES						###
################################################################
INCLUDE(${SOURCE_ROOT}/basics/utilities/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/memory/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/basics/parallel/CMakePackage.txt)
   
#################################################################
###   OWN DEFINES 						###
#################################################################
SET(ALL_SOURCES ${ALL_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp )

SOURCE_GROUP(main FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp )
   
ADD_CXX_FLAGS("/wd4996") #deprecated strcpy...

#################################################################
###   MPI                                                     ###
#################################################################
SET_MPI_STUFF(CAB_MACHINE)
SET(EXECUTABLE_NAME ${EXECUTABLE_NAME}_mpi)

#################################################################
###   EXCECUTABLE						###
#################################################################
ADD_EXECUTABLE(${EXECUTABLE_NAME} ${ALL_SOURCES} )

#################################################################
###   ADDITIONAL LINK LIBRARIES                               ###
#################################################################
IF(ADDITIONAL_LINK_LIBRARIES)
 TARGET_LINK_LIBRARIES(${EXECUTABLE_NAME} ${ADDITIONAL_LINK_LIBRARIES}) 
ENDIF(ADDITIONAL_LINK_LIBRARIES)

#################################################################
###   ADDITIONAL LINK PROPERTIES                              ###
#################################################################
IF(ADDITIONAL_LINK_PROPS)
 SET_TARGET_PROPERTIES(${EXECUTABLE_NAME} PROPERTIES LINK_FLAGS ${ADDITIONAL_LINK_PROPS})
ENDIF(ADDITIONAL_LINK_PROPS)
IF(ADDITIONAL_LINK_PROPS_DEBUG)
 SET_TARGET_PROPERTIES(${EXECUTABLE_NAME} PROPERTIES LINK_FLAGS_DEBUG ${ADDITIONAL_LINK_PROPS_DEBUG})
ENDIF(ADDITIONAL_LINK_PROPS_DEBUG)
IF(ADDITIONAL_LINK_PROPS_RELEASE)
 SET_TARGET_PROPERTIES(${EXECUTABLE_NAME} PROPERTIES LINK_FLAGS_RELEASE ${ADDITIONAL_LINK_PROPS_RELEASE})
ENDIF(ADDITIONAL_LINK_PROPS_RELEASE)



