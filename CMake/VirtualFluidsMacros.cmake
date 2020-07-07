

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
    set( multiValueArgs BUILDTYPE DEPENDS FILES FOLDER EXCLUDE)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    #message("Files: ${ARG_FOLDER}")

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
        file ( GLOB_RECURSE all_files ${VIRTUAL_FLUIDS_GLOB_FILES} )
        set(sourceFiles ${sourceFiles} ${all_files})
    endif()

    if (ARG_EXCLUDE)
        foreach(file_path ${sourceFiles})
            foreach(file_exclude ${ARG_EXCLUDE})
                get_filename_component(file_name ${file_path} NAME)
                if (NOT ${file_name} STREQUAL ${file_exclude})
                    set(new_files ${new_files} ${file_path})
                endif()

            endforeach()
        endforeach()
        set(sourceFiles ${new_files})
    endif()


    includeProductionFiles (${library_name} "${sourceFiles}")

    foreach(X IN LISTS MY_SRCS)
        #message(STATUS "${X}")
    endforeach()


    # SET SOURCE GROUP
    # IF SOURCE GROUP ENABLED
    foreach(source ${sourceFiles})
        #  get_filename_component(source_dir ${source} DIRECTORY)

        if (source_dir)
            #setSourceGroupForFilesIn(${source_dir} ${library_name})
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
        ADD_EXECUTABLE(${library_name} ${MY_SRCS} )
    ELSEIF(${ARG_BUILDTYPE} MATCHES shared)
        ADD_LIBRARY(${library_name} SHARED ${MY_SRCS} )
    ELSEIF(${ARG_BUILDTYPE} MATCHES static)
        ADD_LIBRARY(${library_name} STATIC ${MY_SRCS} )
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

    #MESSAGE (COMPILE FLAGS: ${CAB_ADDTIONAL_COMPILER_FLAGS})
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

    if (NOT ${ARG_BUILDTYPE} MATCHES binary)
      generateExportHeader (${library_name})
    endif()

    target_include_directories(${library_name} PRIVATE ${CMAKE_BINARY_DIR})
    target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${library_name} PRIVATE ${CMAKE_SOURCE_DIR}/src)


endfunction(vf_add_library)

