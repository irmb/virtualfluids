
###############################################################################################################
## Flags ruecksetzen
###############################################################################################################
SET(CAB_ADDTIONAL_COMPILER_FLAGS)

SET(CAB_ADDITIONAL_LINK_FLAGS)
SET(CAB_ADDITIONAL_LINK_FLAGS_DEBUG)
SET(CAB_ADDITIONAL_LINK_FLAGS_RELEASE)

SET(CAB_ADDITIONAL_LINK_LIBRARIES)


LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DSOURCE_ROOT=${VF_ROOT_DIR} )

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


INCLUDE(${VF_CMAKE_DIR}/CMakeSetCompilerFlags.cmake)
INCLUDE(${VF_CMAKE_DIR}/CMakeCompilerMacros.cmake)

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


#################################################################
# ALLGEMEINGUELTIGER STUFF
# CAB_COMPILER setzen und machinespecific configfile laden
###############################################################
SET_CAB_COMPILER()
CHECK_FOR_VARIABLE(CAB_MACHINE "machine name, e.g. ALTIX, ARWEN")
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DCAB_MACHINE_${CAB_MACHINE})
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DCAB_MACHINE=${CAB_MACHINE})
SET(CMAKE_CONFIG_FILE "${VF_CMAKE_DIR}/cmake_config_files/${CAB_MACHINE}.config.cmake")

IF(NOT EXISTS ${CMAKE_CONFIG_FILE})
    status("No configuration file found.")
ELSE()
    status("Load configuration file ${CAB_MACHINE}.config.cmake")
    INCLUDE(${CMAKE_CONFIG_FILE})
ENDIF()
