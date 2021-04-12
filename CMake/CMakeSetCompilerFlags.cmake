
#########################################################################################
## Access the hostname and loads a optional machine file hostname.cmake
#########################################################################################
macro(loadMachineFile)

    site_name(MACHINE_NAME)
    string(TOUPPER  "${MACHINE_NAME}" MACHINE_NAME)

    set(BUILD_MACHINE_FILE_PATH "${VF_CMAKE_DIR}/cmake_config_files")

    set(MACHINE_FILE "${BUILD_MACHINE_FILE_PATH}/${MACHINE_NAME}.config.cmake")

    IF(NOT EXISTS ${MACHINE_FILE})
        status("No configuration file found: ${MACHINE_FILE}.")
    ELSE()
        status("Load configuration file: ${MACHINE_FILE}")
        include(${MACHINE_FILE})
    ENDIF()

endmacro()


################################################################
###               SET_COMPILER_SPECIFIC_FLAGS                ###
###  determines compiler flags variables                     ###
################################################################
macro(loadCompilerFlags)

  SET(CS_COMPILER_FLAGS_CXX "")
  SET(CS_COMPILER_FLAGS_CXX_DEBUG "")
  SET(CS_COMPILER_FLAGS_CXX_RELEASE "")

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
function(addAdditionalFlags library_name)

    status_lib("additional compiler flags CXX: ${CS_COMPILER_FLAGS_CXX}")
    status_lib("additional compiler flags CXX debug: ${CS_COMPILER_FLAGS_CXX_DEBUG}")
    status_lib("additional compiler flags CXX release: ${CS_COMPILER_FLAGS_CXX_RELEASE}")
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
    foreach(flag IN LISTS CS_COMPILER_FLAGS_CXX)
        target_compile_options(${library_name} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${flag}>")
        if(MSVC)
            target_compile_options(${library_name} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${flag}>")
        endif()
    endforeach()

    foreach(flag IN LISTS CS_COMPILER_FLAGS_CXX_DEBUG)
        target_compile_options(${library_name} PRIVATE "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:${flag}>")
    endforeach()

    foreach(flag IN LISTS CS_COMPILER_FLAGS_CXX_RELEASE)
        target_compile_options(${library_name} PRIVATE "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:${flag}>")
    endforeach()

endfunction()