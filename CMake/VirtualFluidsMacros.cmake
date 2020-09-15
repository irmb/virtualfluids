#################################################################################
#   _    ___      __              __________      _     __
# | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
# | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
# | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
# |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
#
#################################################################################

function(status msg)
    message(STATUS "  VF - ${msg}")
endfunction()

#################################################################################
## set global project file endings
#################################################################################
set (VIRTUAL_FLUIDS_GLOB_FILES
        *.cpp
        *.c
        *.h
        *.cu
        *.cuh
        *.hpp
        CACHE INTERNAL "File endings to glob for source files" )


#################################################################################
## Sets the library name to the current folder name.
## output parameter: library_name
#################################################################################
function (vf_get_library_name library_name)
    get_filename_component(library_name_out ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    set(${library_name} ${library_name_out} PARENT_SCOPE)
endfunction()


#################################################################################
## Sets the library test name to the current folder name + Tests.
## output parameter: library_test_name
#################################################################################
function (vf_get_library_test_name library_test_name)
    vf_get_library_name (folder_name)
    set (${library_test_name} ${library_name}Tests PARENT_SCOPE)
endfunction()


#################################################################################
## Add a target, link the libraries and add the compiler flags to the target
##
## parameter:
## NAME      - Name of the target. If not passed the name is vf_get_library_name().
## BUILDTYPE - STATIC; SHARED; EXECUTABLE
## DEPENDS   - libraries to link
## FILES     - adds these files to the target
## FOLDER    - adds all files in these folders to the targets
## EXCLUDE   - excludes these files from the target
##
## note: If no files and folders are passed, all files from the level of current
##       CMakeLists.txt are added recursively.
##
#################################################################################
function(vf_add_library)

    set( options )
    set( oneValueArgs )
    set( multiValueArgs NAME BUILDTYPE DEPENDS FILES FOLDER EXCLUDE)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if(${ARG_NAME})
        set(library_name ${ARG_NAME})
    else()
        vf_get_library_name (library_name)
    endif()
    status("Configuring the target: ${library_name} (type=${ARG_BUILDTYPE})...")


    collectFiles(sourceFiles "${ARG_FILES}" "${ARG_FOLDER}" "${ARG_EXCLUDE}")

    includeProductionFiles (${library_name} "${sourceFiles}")



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
        groupTarget(${library_name} ${appFolder})
    ELSEIF(${ARG_BUILDTYPE} MATCHES shared)
        ADD_LIBRARY(${library_name} SHARED ${MY_SRCS} )
        groupTarget(${library_name} ${libraryFolder})
    ELSEIF(${ARG_BUILDTYPE} MATCHES static)
        ADD_LIBRARY(${library_name} STATIC ${MY_SRCS} )
        groupTarget(${library_name} ${libraryFolder})
    ELSE()
        MESSAGE(FATAL_ERROR "build_type=${ARG_BUILDTYPE} doesn't match BINARY, SHARED or STATIC")
    ENDIF()

    # Set the output directory for build artifacts
    set_target_properties(${library_name}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

    #################################################################
    ###   ADDITIONAL LINK LIBRARIES                               ###
    #################################################################
    status("Link Depending Libraries: ${ARG_DEPENDS}")
    if (ARG_DEPENDS)
        target_link_libraries(${library_name} PRIVATE ${ARG_DEPENDS})
    endif()

    #################################################################
    ###   COMPILER Flags                                          ###
    #################################################################
    ADD_COMPILER_FLAGS_TO_PROJECT(${library_name})
    status("compiler flags for compiler ${CAB_COMPILER} on machine ${CAB_MACHINE} for project ${project_name} (${ARG_BUILDTYPE}) have been configured")
    status("compiler flags: ${CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS}")
    status("additional compiler flags: ${CAB_ADDTIONAL_COMPILER_FLAGS}")

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
    status("additional linker flags: ${CAB_ADDITIONAL_LINK_FLAGS}")

    IF(CAB_ADDITIONAL_LINK_FLAGS)
        ADD_TARGET_PROPERTIES(${library_name} LINK_FLAGS ${CAB_ADDITIONAL_LINK_FLAGS})
    ENDIF()
    IF(CAB_ADDITIONAL_LINK_FLAGS_DEBUG)
        ADD_TARGET_PROPERTIES(${library_name} LINK_FLAGS_DEBUG ${CAB_ADDITIONAL_LINK_FLAGS_DEBUG})
    ENDIF()
    IF(CAB_ADDITIONAL_LINK_FLAGS_RELEASE)
        ADD_TARGET_PROPERTIES(${library_name} LINK_FLAGS_RELEASE ${CAB_ADDITIONAL_LINK_FLAGS_RELEASE})
    ENDIF()


    if (NOT ${ARG_BUILDTYPE} MATCHES binary)
      generateExportHeader (${library_name})
    endif()

    target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${library_name} PRIVATE ${CMAKE_BINARY_DIR})
    target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR})
    target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/gpu)
    target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/cpu)


    #status("... configuring target: ${library_name} (type=${ARG_BUILDTYPE}) done")

endfunction()



#################################################################################
## Add a test executable corresponding to the added target.
## Must be called after vf_add_library().
## The name of the test executable is: vf_get_library_name()Tests
##
## Precondition: BUILD_VF_UNIT_TESTS needs to be ON
#################################################################################
function(vf_add_tests)

    if (NOT BUILD_VF_UNIT_TESTS)
        return()
    endif()

    # get the test library name
    vf_get_library_test_name(library_test_name)
    vf_get_library_name (folder_name)

    status("Configuring test executable: ${library_test_name}")

    # set test files to MY_SRCS
    file ( GLOB_RECURSE all_files ${VIRTUAL_FLUIDS_GLOB_FILES} )
    includeTestFiles (${folder_name} "${all_files}")

    # add the target
    add_executable(${library_test_name} ${MY_SRCS})
    groupTarget (${library_test_name} ${testFolder})

    # Set the output directory for build artifacts
    set_target_properties(${library_test_name}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

    # link tested library
    target_link_libraries(${library_test_name} PRIVATE ${folder_name})

    # link tested library
    target_include_directories(${library_test_name} PRIVATE ${CMAKE_BINARY_DIR})
    target_include_directories(${library_test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${library_test_name} PRIVATE ${VF_SRC_DIR})

    # link googlemock
    linkGMOCK()

endfunction()

#################################################################################
## group target for VisualStudio
#################################################################################
function(groupTarget targetName folderName)
    set_property( TARGET  ${targetName}  PROPERTY  FOLDER  ${folderName} )
endfunction(groupTarget)