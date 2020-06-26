

set (VIRTUAL_FLUIDS_GLOB_FILES
        *.cpp
        *.c
        *.h
        *.cu
        CACHE INTERNAL "File endings to glob for source files" )



function (vf_get_library_name library_name)
    get_filename_component(library_name_out ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    set(${library_name} ${library_name_out} PARENT_SCOPE)
endfunction(vf_get_library_name)



function(vf_add_library)
    message("Start new Cmake")

    set( options )
    set( oneValueArgs )
    set( multiValueArgs DEPENDS FOLDER)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    #message("${ARG_FOLDER}")

    vf_get_library_name (library_name)
    message("${library_name}")

    if (ARG_FOLDER)
        foreach(folder ${ARG_FOLDER})
            foreach(file ${VIRTUAL_FLUIDS_GLOB_FILES})
                set (filePath ${folder}/${file})
                #message("${filePath}")
                file ( GLOB part_files ${filePath} )
                set(sourceFiles ${sourceFiles} ${part_files})
            endforeach()
        endforeach()

        else ()
     file ( GLOB_RECURSE sourceFiles ${VIRTUAL_FLUIDS_GLOB_FILES} )
endif()

    foreach(X IN LISTS sourceFiles)
       #message(STATUS "${X}")
    endforeach()


    set (build_type STATIC)


  MESSAGE(STATUS "configuring ${library_name} (type=${build_type})...")


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
  #SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")

  #################################################################
  ###   EXCECUTABLE                                             ###
  #################################################################
  IF(${build_type} MATCHES BINARY)
      ADD_EXECUTABLE(${library_name} ${sourceFiles} )
  ELSEIF(${build_type} MATCHES SHARED)
      ADD_LIBRARY(${library_name} SHARED ${sourceFiles} )
  ELSEIF(${build_type} MATCHES STATIC)
      ADD_LIBRARY(${library_name} STATIC ${sourceFiles} )
  ELSE()
      MESSAGE(FATAL_ERROR "build_type=${build_type} doesn't match BINARY, SHARED or STATIC")
  ENDIF()

  #################################################################
  ###   ADDITIONAL LINK LIBRARIES                               ###
  #################################################################
    message("Link Depending Libraries: ${ARG_DEPENDS}")
    if (ARG_DEPENDS)
    TARGET_LINK_LIBRARIES(${library_name} ${ARG_DEPENDS})
    endif()

  IF(CAB_ADDITIONAL_LINK_LIBRARIES)
      TARGET_LINK_LIBRARIES(${library_name} ${CAB_ADDITIONAL_LINK_LIBRARIES})
  ENDIF()

  #################################################################
  ###   COMPILER Flags                                          ###
  #################################################################
    message (${CAB_COMPILER})
  ADD_COMPILER_FLAGS_TO_PROJECT(${CAB_COMPILER}  ${library_name} "CXX" ${build_type})
  MESSAGE(STATUS "compiler flags for compiler ${CAB_COMPILER} on machine ${CAB_MACHINE} for project ${project_name} (${build_type}) have been configured")

  IF(CAB_ADDTIONAL_COMPILER_FLAGS)
      ADD_TARGET_PROPERTIES(${library_name} COMPILE_FLAGS ${CAB_ADDTIONAL_COMPILER_FLAGS})
  ENDIF()
  IF(CAB_ADDTIONAL_COMPILER_FLAGS_DEBUG)
      MESSAGE(FATAL_ERROR "COMPILE_FLAGS_DEBUG_<CONFIG> not supported by cmake yet :-(")
      ADD_TARGET_PROPERTIES(${library_name} COMPILE_FLAGS_DEBUG ${CAB_ADDTIONAL_COMPILER_FLAGS_DEBUG})
  ENDIF()
  IF(CAB_ADDTIONAL_COMPILER_FLAGS_RELEASE)
      MESSAGE(FATAL_ERROR "COMPILE_FLAGS_<CONFIG> not supported by cmake yet :-(")
      ADD_TARGET_PROPERTIES(${library_name} COMPILE_FLAGS_RELEASE ${CAB_ADDTIONAL_COMPILER_FLAGS_RELEASE})
  ENDIF()

  #################################################################
  ###   ADDITIONAL LINK PROPERTIES                              ###
  #################################################################
  IF(CAB_ADDITIONAL_LINK_FLAGS)
      ADD_TARGET_PROPERTIES(${library_name} LINK_FLAGS ${CAB_ADDITIONAL_LINK_FLAGS})
  ENDIF()
  IF(CAB_ADDITIONAL_LINK_FLAGS_DEBUG)
      ADD_TARGET_PROPERTIES(${library_name} LINK_FLAGS_DEBUG ${CAB_ADDITIONAL_LINK_FLAGS_DEBUG})
  ENDIF()
  IF(CAB_ADDITIONAL_LINK_FLAGS_RELEASE)
      ADD_TARGET_PROPERTIES(${library_name} LINK_FLAGS_RELEASE ${CAB_ADDITIONAL_LINK_FLAGS_RELEASE})
  ENDIF()

  SET(project_name ${library_name} CACHE STRING "name of binary")

  MESSAGE(STATUS "configuring ${library_name} (type=${build_type})... done")


endfunction(vf_add_library)