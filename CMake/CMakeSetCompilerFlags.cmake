################################################################
###               SET_COMPILER_SPECIFIC_FLAGS                ###
###  determines compiler flags variables                     ###
################################################################
MACRO(LOAD_COMPILER_FLAGS_FROM_FILE)

  SET(CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "")
  SET(CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "")
  SET(CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "")

   # https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_<LANG>_COMPILER_ID

   IF( SPECIFIC_COMPILER_FLAG_FILE )
       include( ${SPECIFIC_COMPILER_FLAG_FILE})
   ELSEIF( EXISTS "${VF_CMAKE_DIR}/compilerflags/${CMAKE_CXX_COMPILER_ID}.cmake" )
       status("Load compiler file: ${CMAKE_CXX_COMPILER_ID}.cmake")
	   include(${VF_CMAKE_DIR}/compilerflags/${CMAKE_CXX_COMPILER_ID}.cmake)
	ELSE()
	   MESSAGE(FATAL_ERROR "compiler=${CMAKE_CXX_COMPILER_ID} seems to be a not supported compiler")
	ENDIF()

ENDMACRO()

################################################################
###             ADD_COMPILER_FLAGS_TO_PROJECT                ###
################################################################
MACRO(ADD_COMPILER_FLAGS_TO_PROJECT project_name)

   LOAD_COMPILER_FLAGS_FROM_FILE()

   #workaround fuer itanium processoren
   IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")
      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS -D_M_IA64)
      LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS   -D_M_IA64)
   ENDIF()


   foreach(flag IN LISTS CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS)
       target_compile_options(${project_name} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${flag}>")
   endforeach()

   foreach(flag IN LISTS CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG)
       target_compile_options(${project_name} PRIVATE "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:${flag}>")
   endforeach()

   foreach(flag IN LISTS CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE)
       target_compile_options(${project_name} PRIVATE "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:${flag}>")
   endforeach()


ENDMACRO(ADD_COMPILER_FLAGS_TO_PROJECT project_name)