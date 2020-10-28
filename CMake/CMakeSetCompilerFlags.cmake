
###############################################################
# set hostname -> CAB_MACHINE and load an optional config file
###############################################################
macro(loadMachineFile)

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

    LIST(APPEND VF_COMPILER_DEFINITION CAB_MACHINE=${CAB_MACHINE})
    SET(CMAKE_CONFIG_FILE "${VF_CMAKE_DIR}/cmake_config_files/${CAB_MACHINE}.config.cmake")

    IF(NOT EXISTS ${CMAKE_CONFIG_FILE})
        status("No configuration file found for machine: ${CAB_MACHINE}.")
    ELSE()
        status("Load configuration file ${CAB_MACHINE}.config.cmake")
        include(${CMAKE_CONFIG_FILE})
    ENDIF()

endmacro()


################################################################
###               SET_COMPILER_SPECIFIC_FLAGS                ###
###  determines compiler flags variables                     ###
################################################################
macro(loadCompilerFlags)

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

endmacro()

################################################################
###             ADD_COMPILER_FLAGS_TO_PROJECT                ###
################################################################
function(addAdditionalFlags project_name)

    status_lib("additional compiler flags CXX: ${CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS}")
    status_lib("additional compiler flags CXX debug: ${CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG}")
    status_lib("additional compiler flags CXX release: ${CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE}")
    status_lib("additional compiler definitions: ${VF_COMPILER_DEFINITION}")
    status_lib("additional linker flags: ${VF_LINK_OPTIONS}")

    # compile definitions
    foreach(flag IN LISTS VF_COMPILER_DEFINITION)
        target_compile_definitions(${library_name} PRIVATE ${flag})
    endforeach()

    # link options
    foreach(flag IN LISTS VF_LINK_OPTIONS) #TODO: check what happens when lib is static
        target_link_options(${library_name} PRIVATE ${flag})
    endforeach()

    # compile options
    foreach(flag IN LISTS CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS)
        target_compile_options(${project_name} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${flag}>")
    endforeach()

    foreach(flag IN LISTS CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG)
        target_compile_options(${project_name} PRIVATE "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:${flag}>")
    endforeach()

    foreach(flag IN LISTS CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE)
        target_compile_options(${project_name} PRIVATE "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:${flag}>")
    endforeach()

endfunction()