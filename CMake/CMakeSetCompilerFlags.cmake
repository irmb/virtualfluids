################################################################
###               SET_COMPILER_SPECIFIC_FLAGS                ###
###  determines compiler flags variabels                     ###
###  compiler_type: e.g. msvc9_x64                           ###
###  build_type  :    BINARY, STATIC, SHARED                 ###
################################################################
MACRO(SET_COMPILER_SPECIFIC_FLAGS compiler_type build_type)
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

   ###############################################################################################################
   ## ggf. spezielles compiler flag file lesen
   ###############################################################################################################
   IF( SPECIFIC_COMPILER_FLAG_FILE )
      INCLUDE( ${SPECIFIC_COMPILER_FLAG_FILE})
   ###############################################################################################################
   ## standard compiler flags
   ###############################################################################################################
   ELSEIF( EXISTS "${VF_CMAKE_DIR}/compilerflags/${CAB_COMPILER}.cmake" )
       status("Load compiler file: ${CAB_COMPILER}.cmake")
	   INCLUDE( ${VF_CMAKE_DIR}/compilerflags/${CAB_COMPILER}.cmake)
	###############################################################################################################
	## unknown compiler
	###############################################################################################################
	ELSE()
	   #MESSAGE(FATAL_ERROR "CAB_COMPILER=${CAB_COMPILER} seems to be a not supported compiler")
	   #MESSAGE(WARNING "CAB_COMPILER=${CAB_COMPILER} seems to be a not supported compiler; set to generic")
	   SET(CAB_COMPILER "gccGeneric")
	   INCLUDE( ${VF_CMAKE_DIR}/compilerflags/${CAB_COMPILER}.cmake)
	ENDIF()
   

   ###############################################################################################################
	#64 Bit compilation??
   ###############################################################################################################
   IF(NOT DEFINED USE_64BIT_COMPILER_OPTIONS) 
      IF(MSVC AND NOT CMAKE_CL_64)      
        SET(OPTION64 OFF)
      ELSEIF(MSVC AND CMAKE_CL_64)      
        SET(OPTION64 ON)
      ELSE()
         IS_64BIT_SYSTEM( IS64BITSYSTEM )
         IF(IS64BITSYSTEM STREQUAL "TRUE")
            SET(OPTION64 ON)
         ELSE()
            SET(OPTION64 OFF)
         ENDIF()
      ENDIF()
   ENDIF()

   OPTION(USE_64BIT_COMPILER_OPTIONS "set 64 bit compiler flags"  ${OPTION64})

   ###############################################################################################################
	# set flags
   ###############################################################################################################
    IF(USE_64BIT_COMPILER_OPTIONS)
          SET_COMPILER_SPECIFIC_FLAGS_INTERN( ${build_type} 1)
    else()
          SET_COMPILER_SPECIFIC_FLAGS_INTERN( ${build_type} 0)
    endif()
  
ENDMACRO(SET_COMPILER_SPECIFIC_FLAGS compiler_type build_type)

################################################################
###             ADD_COMPILER_FLAGS_TO_PROJECT                ###
###  adds COMPILER_FLGAS TO project                          ###
###  project_language:    CXX , C                            ###
################################################################
MACRO(ADD_COMPILER_FLAGS_TO_PROJECT compiler_type project_name project_language build_type)

   IF(NOT ${project_language} MATCHES "C")
      IF(NOT ${project_language} MATCHES "CXX")
         MESSAGE(FATAL_ERROR "project_language must be CXX or C")
      ENDIF()
   ENDIF()

   ################################################################
   # SET_COMPILER_SPECIFIC_FLAGS
   ################################################################
   SET_COMPILER_SPECIFIC_FLAGS( ${compiler_type} ${build_type} )

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


ENDMACRO(ADD_COMPILER_FLAGS_TO_PROJECT compiler_type project_name project_language build_type)