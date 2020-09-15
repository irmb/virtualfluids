################################################################
###               SET_COMPILER_SPECIFIC_FLAGS                ###
###  determines compiler flags variables                     ###
################################################################
MACRO(SET_COMPILER_SPECIFIC_FLAGS)
   IF(NOT CMAKE_CXX_COMPILER)
      MESSAGE(FATAL_ERROR "before SET_CAB_COMPILER-Macro PROJECT-Macro has to be called")
   ENDIF()

  ###############################################################################################################
  ## Flags ruecksetzen
  ###############################################################################################################
  SET(CAB_COMPILER_ADDITIONAL_LINK_PROPS "")
  SET(CAB_COMPILER_ADDITIONAL_LINK_PROPS_DEBUG "")
  SET(CAB_COMPILER_ADDITIONAL_LINK_PROPS_RELEASE "")
  
  SET(CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "")
  SET(CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "")
  SET(CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "")

  SET(CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "")
  SET(CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS_DEBUG "")
  SET(CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS_RELEASE "")


   # https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_<LANG>_COMPILER_ID

   IF( SPECIFIC_COMPILER_FLAG_FILE )
       include( ${SPECIFIC_COMPILER_FLAG_FILE})
   ELSEIF( EXISTS "${VF_CMAKE_DIR}/compilerflags/${CMAKE_CXX_COMPILER_ID}.cmake" )
       status("Load compiler file: ${CMAKE_CXX_COMPILER_ID}.cmake")
	   include(${VF_CMAKE_DIR}/compilerflags/${CMAKE_CXX_COMPILER_ID}.cmake)
	ELSE()
	   MESSAGE(FATAL_ERROR "compiler=${CMAKE_CXX_COMPILER_ID} seems to be a not supported compiler")
	ENDIF()

ENDMACRO(SET_COMPILER_SPECIFIC_FLAGS)

################################################################
###             ADD_COMPILER_FLAGS_TO_PROJECT                ###
###  adds COMPILER_FLGAS TO project                          ###
################################################################
MACRO(ADD_COMPILER_FLAGS_TO_PROJECT project_name)

   ################################################################
   # SET_COMPILER_SPECIFIC_FLAGS
   ################################################################
   SET_COMPILER_SPECIFIC_FLAGS()

   #workaround fuer itanium processoren
   IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")
      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS -D_M_IA64)
      LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS   -D_M_IA64)
   ENDIF()

   ################################################################
   # LINKER PROPS
   ################################################################
   status("additional linker probs: ${CAB_COMPILER_ADDITIONAL_LINK_PROPS}")

   IF(CAB_COMPILER_ADDITIONAL_LINK_PROPS)
     ADD_TARGET_PROPERTIES(${project_name} LINK_FLAGS ${CAB_COMPILER_ADDITIONAL_LINK_PROPS}) 
   ENDIF()
   IF(CAB_COMPILER_ADDITIONAL_LINK_PROPS_DEBUG)
     ADD_TARGET_PROPERTIES(${project_name} LINK_FLAGS ${CAB_COMPILER_ADDITIONAL_LINK_PROPS_DEBUG}) 
   ENDIF()
   IF(CAB_COMPILER_ADDITIONAL_LINK_PROPS_RELEASE)
     ADD_TARGET_PROPERTIES(${project_name} LINK_FLAGS ${CAB_COMPILER_ADDITIONAL_LINK_PROPS_RELEASE})
   ENDIF()

   ################################################################
   # COMPILER FLAGS
   ################################################################
   #message (COMPILE FLAGS INTERN: ${CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS})

    # TODO: Clean this up!!
   IF(CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS)
       #message (COMPILE FLAGS INTERN: ${CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS})
       foreach(flag IN LISTS CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS)
           #message(compiler option: ${flag})
           target_compile_options(${project_name} PRIVATE "$<$<COMPILE_LANGUAGE:${project_language}>:${flag}>")
       endforeach()
       #get_target_property(var ${project_name} COMPILE_OPTIONS)
       #message(set compile options: ${var})

       #add_custom_command(TARGET ${project_name} POST_BUILD COMMAND echo built with the flags: ${var})
       #ADD_TARGET_PROPERTIES(${project_name} COMPILE_FLAGS ${CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS})
       #target_compile_options (${project_name} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS}>)
   ENDIF()
   IF(CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS_DEBUG)
     MESSAGE(STATUS "ADD_COMPILER_FLAGS_TO_PROJECT: sorry, a long as CMake has no support for COMPILE_FLAGS_<CONFIG> -> DEBUG flags are neglected")
     #ADD_TARGET_PROPERTIES(${project_name} COMPILE_FLAGS_DEBUG ${CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS_DEBUG})
   ENDIF()
   IF(CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS_RELEASE)
     MESSAGE(STATUS "ADD_COMPILER_FLAGS_TO_PROJECT: sorry, a long as CMake has no support for COMPILE_FLAGS_<CONFIG> -> RELEASE flags are set for RELEASE AND DEBUG")
     ADD_TARGET_PROPERTIES(${project_name} COMPILE_FLAGS ${CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS_RELEASE}) 
     #ADD_TARGET_PROPERTIES(${project_name} COMPILE_FLAGS_RELEASE ${CAB_COMPILER_ADDTIONAL_${project_language}_COMPILER_FLAGS_RELEASE})
   ENDIF()


ENDMACRO(ADD_COMPILER_FLAGS_TO_PROJECT project_name)