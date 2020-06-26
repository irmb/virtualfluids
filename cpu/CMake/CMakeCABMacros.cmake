###############################################################
# Aktivieren von IF(ARG)...ELSE()...ENDIF() in CMake
###############################################################
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS 1)

###############################################################################################################
## Flags ruecksetzen
###############################################################################################################
SET(CAB_ADDTIONAL_COMPILER_FLAGS)
#because of the fact that CMake does not support COMPILER_FLAGS_<CONFIG> right know we cannot use these options
#SET(CAB_ADDTIONAL_COMPILER_FLAGS_DEBUG)
#SET(CAB_ADDTIONAL_COMPILER_FLAGS_RELEASE)

SET(CAB_ADDITIONAL_LINK_FLAGS)
SET(CAB_ADDITIONAL_LINK_FLAGS_DEBUG)
SET(CAB_ADDITIONAL_LINK_FLAGS_RELEASE)

SET(CAB_ADDITIONAL_LINK_LIBRARIES)


###############################################################
# SOURCE_DIR variable wird beim einbinden automatisch gesetzt - ebenso das include dir!!!
# SOURCE_DIR wird dem Projekt als Standardincludepfad hinzugefuegt
###############################################################
#CMakeCABMacros.txt liegt direkt im source-Ordner -> pfad == SOURCE_ROOT
GET_FILENAME_COMPONENT( SOURCE_ROOT  ${CMAKE_CURRENT_LIST_FILE} PATH) 
STRING(REGEX REPLACE "(.*)/CMake" "\\1" SOURCE_ROOT "${SOURCE_ROOT}" )
INCLUDE_DIRECTORIES(${SOURCE_ROOT})
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DSOURCE_ROOT=${SOURCE_ROOT} )

###############################################################
# hostename ermitteln -> CAB_MACHINE
###############################################################
IF(NOT CAB_MACHINE)
   SET(CAB_MACHINE $ENV{CAB_MACHINE})

   IF( CAB_MACHINE )
     STRING(TOUPPER  "${CAB_MACHINE}" CAB_MACHINE)
   ELSE()
     EXECUTE_PROCESS( COMMAND hostname OUTPUT_VARIABLE CAB_MACHINE)
     STRING(REGEX REPLACE "[ ]*([A-Za-z0-9]+).*[\\\\n]*" "\\1" CAB_MACHINE "${CAB_MACHINE}" )
     STRING(TOUPPER  "${CAB_MACHINE}" CAB_MACHINE)
   ENDIF()
ENDIF()


###############################################################
# WITH_SUBFOLDERS_FOR_SG erstellen, wenn noch nicht vorhanden
# ist diese auf true wird bei SOURCE_GROUP in VS unterordner erzeugt
# dies funzt erst ab CMake-2.5 vernuenftig
###############################################################
IF(NOT WITH_SUBFOLDERS_FOR_SG)
   SET(WITH_SUBFOLDERS_FOR_SG FALSE)
ENDIF()


############################################################################
############################################################################
##                  M      A      C      R      O      S                  ##
##                  M      A      C      R      O      S                  ##
##                  M      A      C      R      O      S                  ##
############################################################################
############################################################################
# externe (ACHTUNG: die darin enthaltenen benoetigen teils noch macros die 
# hier im Anschluss folgen
INCLUDE("${SOURCE_ROOT}/CMake/CMakeSetCompilerFlags.cmake")
INCLUDE("${SOURCE_ROOT}/CMake/CMakeCompilerMacros.cmake")

################################################################
###               ADD_TARGET_PROPERTIES                      ###
################################################################
MACRO(ADD_TARGET_PROPERTIES target property)

   SET(property_values ${ARGN})
   
   # vorhandene properties holen
   GET_TARGET_PROPERTY(TEMP_VAR ${target} ${property}) 
   IF(TEMP_VAR)
      LIST(APPEND property_values ${TEMP_VAR})     
   ENDIF() 

   #STRING(REGEX REPLACE ";" " " property_values ${property_values})
   # doppelte Eintraege loeschen
   SEPARATE_ARGUMENTS(property_values)
   LIST(REMOVE_DUPLICATES property_values)

   #aus Liste wieder einen String basteln (geht nich tmit regex replace...)
   SET(new_property_values)
   FOREACH(p ${property_values})
      SET(new_property_values "${p} ${new_property_values}")
   ENDFOREACH()

   #property setzen
   SET_TARGET_PROPERTIES(${target} PROPERTIES ${property} ${new_property_values})

   #GET_TARGET_PROPERTY(TEMP_VAR ${target} ${property}) 
   #MESSAGE("danach ${target} ${property} ${TEMP_VAR}")

ENDMACRO(ADD_TARGET_PROPERTIES target property)


################################################################
###               CHECK_FOR_VARIABLE                         ###
###  checks for a variable (also env-variables)
###  if not found -> error-message!!!
###  always: cache-entry update
################################################################
MACRO(CHECK_FOR_VARIABLE var)
   #check ob evtl enviromentvariable gesetzt
   IF(NOT DEFINED ${var})  #true if ${var} NOT: empty, 0, N, NO, OFF, FALSE, NOTFOUND, or <variable>-NOTFOUND
     SET(${var} $ENV{${var}})
   ENDIF()

   IF(NOT DEFINED ${var})
      SET(${var} "${var}-NOTFOUND" CACHE STRING "${ARGN}" FORCE)
   ENDIF(NOT DEFINED ${var})

   IF(${var} MATCHES ".*-NOTFOUND")
      MESSAGE(FATAL_ERROR "CHECK_FOR_VARIABLE - error - set ${var}")
   ENDIF()
   
   SET(${var} ${${var}} CACHE STRING "${ARGN}" FORCE) 
ENDMACRO(CHECK_FOR_VARIABLE var)


###############################################################
###   CAB_SOURCE_GROUP( sourceGroupName )                   ###
### creates a source group for the given folder and files.  ### 
###############################################################
MACRO(CAB_SOURCE_GROUP sourceGroupName)
  SET(tempSourceGroupName ${sourceGroupName})
  IF(WITH_SUBFOLDERS_FOR_SG)
    STRING(REGEX REPLACE "/" "\\\\" tempSourceGroupName ${tempSourceGroupName})
  ENDIF()
  SOURCE_GROUP(${tempSourceGroupName} FILES ${ARGN})
ENDMACRO(CAB_SOURCE_GROUP)

#################################################################################
###   COLLECT_PACKAGE_DATA( currentDir  sourceGroupName outFiles)             ###
### collects header and cpp file of current dir and add them to "outfiles"    ###
### all files will be put to the SOURCE_GROUP-folder "sourceGroupName"        ###
### and this one will be with subfolders if  WITH_SUBFOLDERS_FOR_SG==YES      ###
#################################################################################
MACRO(COLLECT_PACKAGE_DATA currentDir sourceGroupName outFiles)
  FILE( GLOB _HEADER_FILES ${currentDir}/*.h   )
  FILE( GLOB _CPP_FILES    ${currentDir}/*.cpp )
  FILE( GLOB _CXX_FILES    ${currentDir}/*.cxx )
  FILE( GLOB _HPP_FILES    ${currentDir}/*.hpp )
  FILE( GLOB _C_FILES      ${currentDir}/*.c   )
  
  IF(CAB_PACKAGE_DEFINTIONS)
     SET_SOURCE_FILES_PROPERTIES( ${CPP_FILES} PROPERTIES COMPILE_FLAGS ${CAB_PACKAGE_DEFINTIONS} )
  ENDIF(CAB_PACKAGE_DEFINTIONS)

  CAB_SOURCE_GROUP( ${sourceGroupName} ${_HEADER_FILES} ${_CPP_FILES} ${_CXX_FILES} ${_HPP_FILES} ${_C_FILES} )

  #SET( ${outFiles} ${${outFiles}} ${_HEADER_FILES} ${_CPP_FILES} ${_CXX_FILES} ${_HPP_FILES} ${_C_FILES} ) 
  LIST(APPEND ${outFiles} ${_HEADER_FILES} ${_CPP_FILES} ${_CXX_FILES} ${_HPP_FILES} ${_C_FILES} )
ENDMACRO(COLLECT_PACKAGE_DATA  currentDir sourceGroupName sourceGroupWithSubfolders outFiles)


##################################################################################################################
###   COLLECT_PACKAGE_DATA_WITH_OPTION( currentDir  outFiles [outOption] [outSourceGroupName]) ###
### collects header and cpp file of current dir and add them to "outfiles"                                     ###
### all files will be put to the SOURCE_GROUP-folder "sourceGroupName"                                         ###
### and this one will be with subfolders if  WITH_SUBFOLDERS_FOR_SG==YES                                       ###
##################################################################################################################
MACRO(COLLECT_PACKAGE_DATA_WITH_OPTION currentDir outFiles)
  STRING(REGEX REPLACE "(.*)/source/(.*)" "\\2" SOURCE_GROUP_NAME "${currentDir}")
  STRING(REGEX REPLACE   "/" "_" OPTION_LABEL "${SOURCE_GROUP_NAME}")
  STRING(REGEX REPLACE   ":" "" OPTION_LABEL "${OPTION_LABEL}")
  STRING(TOUPPER ${OPTION_LABEL} OPTION_LABEL)  
      
  SET(OPTION_LABEL "BUILD_${OPTION_LABEL}")
  OPTION(${OPTION_LABEL} "${currentDir}" ON)
  
  IF( ${OPTION_LABEL} ) 
     COLLECT_PACKAGE_DATA( ${currentDir} ${SOURCE_GROUP_NAME} ${outFiles})
  ENDIF(${OPTION_LABEL})

  IF(${ARGC} GREATER 2)
    SET( ${ARGV2} ${OPTION_LABEL} )
  ENDIF()

  IF(${ARGC} GREATER 3)
     SET( ${ARGV3} ${SOURCE_GROUP_NAME} )
  ENDIF()
  
ENDMACRO(COLLECT_PACKAGE_DATA_WITH_OPTION  currentDir outFiles)


#################################################################
###   GET_DIRECTORY_FROM_ENV( var env_var [CACHE] [STRING])   ###
### if enn_var exists the value with corrected slashes will   ###
### be stored in var.					      ###
### if optional CACHE is activated the var will be stored as  ###
### cache variable. optional you can use a status bar string   ###
#################################################################
MACRO(GET_DIRECTORY_FROM_ENV var env_var)
  SET(${var} $ENV{${env_var}})
  IF(${var}) 
    STRING(REGEX REPLACE "\\\\" "/" ${var} ${${var}})   # "\" --> "/" 
    IF(${ARGC} EQUAL 3 AND ${ARGV2} MATCHES "CACHE")
       SET(${var} ${${var}} CACHE PATH "" FORCE)
    ENDIF(${ARGC} EQUAL 3 AND ${ARGV2} MATCHES "CACHE")
    IF(${ARGC} EQUAL 4 AND ${ARGV2} MATCHES "CACHE")
       SET(${var} ${${var}} CACHE PATH "${ARGV3}" FORCE)
    ENDIF(${ARGC} EQUAL 4 AND ${ARGV2} MATCHES "CACHE")
  ENDIF(${var})
ENDMACRO(GET_DIRECTORY_FROM_ENV var env_var)

#################################################################
###   GET_DIRECTORY_FROM_VAR( var  [CACHE] [STRING])          ###
### if optional CACHE is activated the var will be stored as  ###
### cache variable. optional you can use a status bar string   ###
#################################################################
MACRO(GET_DIRECTORY_FROM_VAR var )
  IF(${var}) 
    STRING(REGEX REPLACE "\\\\" "/" ${var} ${${var}})   # "\" --> "/" 
    IF(${ARGC} EQUAL 2 AND ${ARGV1} MATCHES "CACHE")
       SET(${var} ${${var}} CACHE PATH "" FORCE)
    ENDIF(${ARGC} EQUAL 2 AND ${ARGV1} MATCHES "CACHE")
    IF(${ARGC} EQUAL 3 AND ${ARGV1} MATCHES "CACHE")
       SET(${var} ${${var}} CACHE PATH "${ARGV2}" FORCE)
    ENDIF(${ARGC} EQUAL 3 AND ${ARGV1} MATCHES "CACHE")
  ENDIF(${var})
ENDMACRO(GET_DIRECTORY_FROM_VAR var env_var)


#################################################################
###   FINAL MACROS TO GENERATE PROJECT FILES                  ###
###   project_name: name of the project                       ###
###   build_type:   BINARY | SHARED | STATIC                  ###
###   optinal:      prefix
#################################################################
MACRO(CREATE_CAB_PROJECT project_name build_type)

   MESSAGE(STATUS "configuring ${project_name} (type=${build_type})...")


   #################################################################
   ###   OS DEFINES                                              ###
   #################################################################
   IF(WIN32)  
      LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__WIN__)
   ELSEIF(APPLE)
      LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__APPLE__)
   ELSEIF(UNIX)
      LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
   ENDIF()

   #################################################################
   ###   ADDITIONAL_MAKE_CLEAN_FILES                             ###
   #################################################################
   SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")
   
   #################################################################
   ###   EXCECUTABLE                                             ###
   #################################################################
   IF(${build_type} MATCHES BINARY)
      ADD_EXECUTABLE(${project_name} ${ALL_SOURCES} )
   ELSEIF(${build_type} MATCHES SHARED)
      ADD_LIBRARY(${project_name} SHARED ${ALL_SOURCES} )
   ELSEIF(${build_type} MATCHES STATIC)
      ADD_LIBRARY(${project_name} STATIC ${ALL_SOURCES} )
   ELSE()
      MESSAGE(FATAL_ERROR "build_type=${build_type} doesn't match BINARY, SHARED or STATIC")
   ENDIF()

   #################################################################
   ###   ADDITIONAL LINK LIBRARIES                               ###
   #################################################################
   IF(CAB_ADDITIONAL_LINK_LIBRARIES)
      TARGET_LINK_LIBRARIES(${project_name} ${CAB_ADDITIONAL_LINK_LIBRARIES}) 
   ENDIF()

   #################################################################
   ###   COMPILER Flags                                          ###
   #################################################################
   ADD_COMPILER_FLAGS_TO_PROJECT(${CAB_COMPILER}  ${project_name} "CXX" ${build_type})
   MESSAGE(STATUS "compiler flags for compiler ${CAB_COMPILER} on machine ${CAB_MACHINE} for project ${project_name} (${build_type}) have been configured")

   IF(CAB_ADDTIONAL_COMPILER_FLAGS)
     ADD_TARGET_PROPERTIES(${project_name} COMPILE_FLAGS ${CAB_ADDTIONAL_COMPILER_FLAGS}) 
   ENDIF()
   IF(CAB_ADDTIONAL_COMPILER_FLAGS_DEBUG)
     MESSAGE(FATAL_ERROR "COMPILE_FLAGS_DEBUG_<CONFIG> not supported by cmake yet :-(")
     ADD_TARGET_PROPERTIES(${project_name} COMPILE_FLAGS_DEBUG ${CAB_ADDTIONAL_COMPILER_FLAGS_DEBUG}) 
   ENDIF()
   IF(CAB_ADDTIONAL_COMPILER_FLAGS_RELEASE)
     MESSAGE(FATAL_ERROR "COMPILE_FLAGS_<CONFIG> not supported by cmake yet :-(")
     ADD_TARGET_PROPERTIES(${project_name} COMPILE_FLAGS_RELEASE ${CAB_ADDTIONAL_COMPILER_FLAGS_RELEASE}) 
   ENDIF()

   #################################################################
   ###   ADDITIONAL LINK PROPERTIES                              ###
   #################################################################
   IF(CAB_ADDITIONAL_LINK_FLAGS)
     ADD_TARGET_PROPERTIES(${project_name} LINK_FLAGS ${CAB_ADDITIONAL_LINK_FLAGS}) 
   ENDIF()
   IF(CAB_ADDITIONAL_LINK_FLAGS_DEBUG)
     ADD_TARGET_PROPERTIES(${project_name} LINK_FLAGS_DEBUG ${CAB_ADDITIONAL_LINK_FLAGS_DEBUG}) 
   ENDIF()
   IF(CAB_ADDITIONAL_LINK_FLAGS_RELEASE)
     ADD_TARGET_PROPERTIES(${project_name} LINK_FLAGS_RELEASE ${CAB_ADDITIONAL_LINK_FLAGS_RELEASE}) 
   ENDIF()

   SET(project_name ${project_name} CACHE STRING "name of binary")

   MESSAGE(STATUS "configuring ${project_name} (type=${build_type})... done")

ENDMACRO(CREATE_CAB_PROJECT project_name build_type)

#################################################################
# ALLGEMEINGUELTIGER STUFF
# CAB_COMPILER setzen und machinespecific configfile laden
###############################################################
SET_CAB_COMPILER()
CHECK_FOR_VARIABLE(CAB_MACHINE "machine name, e.g. ALTIX, ARWEN")
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DCAB_MACHINE_${CAB_MACHINE})
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DCAB_MACHINE=${CAB_MACHINE})
SET(CMAKE_CONFIG_FILE "${SOURCE_ROOT}/CMake/cmake_config_files/${CAB_MACHINE}.config.cmake")

IF(NOT EXISTS ${CMAKE_CONFIG_FILE})
  MESSAGE(FATAL_ERROR "${CMAKE_CONFIG_FILE} does not exists... maybe false CAB_MACHINE = ${CAB_MACHINE}")
ELSE()
  INCLUDE(${CMAKE_CONFIG_FILE})  
ENDIF()