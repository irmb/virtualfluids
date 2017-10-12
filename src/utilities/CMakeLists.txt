
#############################################################
###                PROJECT                                ###
#############################################################
set(projectName "utilities")
project(${projectName})
MESSAGE(STATUS "Build project \"${projectName}\":")
#############################################################
###		GENERAL MACROS          						  ###
#############################################################
include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/BuildTarget.cmake)
include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/BuildTargetUtilities.cmake)
include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/EndingsToCollect.cmake)
include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/general/FileUtilities.cmake)

#############################################################
###		GENERAL SETTINGS         						  ###
#############################################################
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")
sharedLibs()

#############################################################
###		SUBDIRECTORIES									  ###
#############################################################
add_subdirectory(applications/utilities)
add_subdirectory(applications/utilities.Test)