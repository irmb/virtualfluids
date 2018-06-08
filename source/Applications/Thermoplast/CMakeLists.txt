CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(thermoplast)
IF(${USE_DEM_COUPLING})
	INCLUDE(${APPS_ROOT}/IncludsList.cmake) 
	INCLUDE(${SOURCE_ROOT}/DemCoupling/IncludsList.cmake)

	#################################################################
	###   LOCAL FILES                                             ###
	#################################################################
	FILE(GLOB SPECIFIC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h
							 ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
							 ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp  )
	 
	SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
	SOURCE_GROUP(src FILES ${SPECIFIC_FILES})
	  
	SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} VirtualFluids)
	
	#message("CAB_ADDITIONAL_LINK_LIBRARIES: " ${CAB_ADDITIONAL_LINK_LIBRARIES})

	#################################################################
	###   CREATE PROJECT                                          ###
	#################################################################
	CREATE_CAB_PROJECT(thermoplast BINARY)
ENDIF()