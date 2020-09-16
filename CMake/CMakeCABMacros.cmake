
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
    IF(NOT ${var})  #true if ${var} NOT: empty, 0, N, NO, OFF, FALSE, NOTFOUND, or <variable>-NOTFOUND
        SET(${var} $ENV{${var}})
    ENDIF()

    IF(NOT DEFINED ${var})
        SET(${var} "${var}-NOTFOUND" CACHE STRING "${ARGN}" FORCE)
    ENDIF(NOT DEFINED ${var})

    IF(NOT ${var})
        MESSAGE(FATAL_ERROR "CHECK_FOR_VARIABLE - error - set ${var}")
    ENDIF()

    SET(${var} ${${var}} CACHE STRING "${ARGN}" FORCE)
ENDMACRO(CHECK_FOR_VARIABLE var)



INCLUDE(${VF_CMAKE_DIR}/CMakeSetCompilerFlags.cmake)

###############################################################################################################
# Reset the compiler and linker flags
###############################################################################################################
SET(VF_COMPILER_DEFINITION)
SET(VF_LINK_OPTIONS)

SET(CAB_ADDITIONAL_LINK_LIBRARIES)

LIST(APPEND VF_COMPILER_DEFINITION SOURCE_ROOT=${VF_ROOT_DIR} )

#################################################################
###   OS DEFINES                                              ###
#################################################################
IF(WIN32)
    list(APPEND VF_COMPILER_DEFINITION __WIN__)
ELSEIF(UNIX)
    list(APPEND VF_COMPILER_DEFINITION __unix__)
ENDIF()
IF(APPLE)
    list(APPEND VF_COMPILER_DEFINITION __APPLE__)
endif()



###############################################################
# set hostname -> CAB_MACHINE and load an optional config file
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

CHECK_FOR_VARIABLE(CAB_MACHINE "machine name, e.g. ALTIX, ARWEN")
LIST(APPEND VF_COMPILER_DEFINITION CAB_MACHINE=${CAB_MACHINE})
SET(CMAKE_CONFIG_FILE "${VF_CMAKE_DIR}/cmake_config_files/${CAB_MACHINE}.config.cmake")

IF(NOT EXISTS ${CMAKE_CONFIG_FILE})
    status("No configuration file found.")
ELSE()
    status("Load configuration file ${CAB_MACHINE}.config.cmake")
    INCLUDE(${CMAKE_CONFIG_FILE})
ENDIF()
