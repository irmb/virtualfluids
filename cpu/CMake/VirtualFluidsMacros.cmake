

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
    set( multiValueArgs BUILDTYPE DEPENDS FILES FOLDER)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    message("Files: ${ARG_FILES}")

    vf_get_library_name (library_name)

    if (ARG_FILES)
        set(sourceFiles ${sourceFiles} ${ARG_FILES})
    endif()

    if (ARG_FOLDER)
        foreach(folder ${ARG_FOLDER})
            foreach(file ${VIRTUAL_FLUIDS_GLOB_FILES})
                set (filePath ${folder}/${file})
                #message("${filePath}")
                file (GLOB part_files ${filePath} )
                set(sourceFiles ${sourceFiles} ${part_files})
            endforeach()
        endforeach()
    endif()

    if (NOT ARG_FILES AND NOT ARG_FOLDER)
        file ( GLOB_RECURSE sourceFiles ${VIRTUAL_FLUIDS_GLOB_FILES} )
    endif()

    #    foreach(X IN LISTS sourceFiles)
    #      message(STATUS "${X}")
    #    endforeach()


    # SET SOURCE GROUP
    # IF SOURCE GROUP ENABLED
    foreach(source ${sourceFiles})
        get_filename_component(source_dir ${source} DIRECTORY)

        if (source_dir)
          setSourceGroupForFilesIn(${source_dir} ${library_name})
        endif()
    endforeach()



    MESSAGE(STATUS "configuring ${library_name} (type=${ARG_BUILDTYPE})...")


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
    IF(${ARG_BUILDTYPE} MATCHES binary)
        ADD_EXECUTABLE(${library_name} ${sourceFiles} )
    ELSEIF(${ARG_BUILDTYPE} MATCHES shared)
        ADD_LIBRARY(${library_name} shared ${sourceFiles} )
    ELSEIF(${ARG_BUILDTYPE} MATCHES static)
        ADD_LIBRARY(${library_name} STATIC ${sourceFiles} )
    ELSE()
        MESSAGE(FATAL_ERROR "build_type=${ARG_BUILDTYPE} doesn't match BINARY, SHARED or STATIC")
    ENDIF()

    #################################################################
    ###   ADDITIONAL LINK LIBRARIES                               ###
    #################################################################
    message("Link Depending Libraries: ${ARG_DEPENDS}")
    if (ARG_DEPENDS)
        target_link_libraries(${library_name} PRIVATE ${ARG_DEPENDS})
    endif()

    IF(CAB_ADDITIONAL_LINK_LIBRARIES)
        TARGET_LINK_LIBRARIES(${library_name} PRIVATE ${CAB_ADDITIONAL_LINK_LIBRARIES})
    ENDIF()

    #################################################################
    ###   COMPILER Flags                                          ###
    #################################################################
    ADD_COMPILER_FLAGS_TO_PROJECT(${CAB_COMPILER} ${library_name} "CXX" ${ARG_BUILDTYPE})
    MESSAGE(STATUS "compiler flags for compiler ${CAB_COMPILER} on machine ${CAB_MACHINE} for project ${project_name} (${ARG_BUILDTYPE}) have been configured")

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

    #SET(project_name ${library_name} CACHE STRING "name of binary")

    MESSAGE(STATUS "configuring ${library_name} (type=${ARG_BUILDTYPE})... done")

endfunction(vf_add_library)




macro(setSourceGroupForFilesIn package_dir targetName)
    #input: target_name PACKAGE_SRCS
    buildSourceGroup(${targetName} ${package_dir})

    if(isAllTestSuite)
        source_group(${targetName}\\${SOURCE_GROUP} FILES ${source})
    else()
        source_group(${SOURCE_GROUP} FILES ${source})
    endif()
    #output: -
endmacro(setSourceGroupForFilesIn)




macro(buildSourceGroup targetName path)
    #input: targetName (e.g. lib name, exe name)

    unset(SOURCE_GROUP)
    string(REPLACE "/" ";" folderListFromPath ${path})
    set(findTargetName 0)

    foreach(folder ${folderListFromPath})
        if(findTargetName)
            set(SOURCE_GROUP ${SOURCE_GROUP}\\${folder})
        endif()

        if(${folder} STREQUAL ${targetName})
            SET(findTargetName 1)
        endif()
    endforeach()

    if(NOT SOURCE_GROUP)
        set(SOURCE_GROUP "general")
    endif()
    message("Source group: ${SOURCE_GROUP}")
    #output: SOURCE_GROUP
endmacro(buildSourceGroup)