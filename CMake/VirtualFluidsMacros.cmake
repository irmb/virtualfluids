#=======================================================================================
# ____          ____    __    ______     __________   __      __       __        __
# \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
#  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
#   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
#    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
#     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
#      \    \  |    |   ________________________________________________________________
#       \    \ |    |  |  ______________________________________________________________|
#        \    \|    |  |  |         __          __     __     __     ______      _______
#         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
#          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
#           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
#            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
#
#  This file is part of VirtualFluids. VirtualFluids is free software: you can
#  redistribute it and/or modify it under the terms of the GNU General Public
#  License as published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  SPDX-License-Identifier: GPL-3.0-or-later
#  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
#
#! \author Soeren Peters
#=======================================================================================

function(status msg)
    message(STATUS "  VF: ${msg}")
endfunction()

function(status_lib msg)
    message(DEBUG "    ${msg}")
endfunction()

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
    set (${library_test_name} ${folder_name}Tests PARENT_SCOPE)
endfunction()

#################################################################################
## Add a target, link the libraries and add the compiler flags to the target
##
## parameter:
## NAME      - Name of the target. If not passed the name is vf_get_library_name().
## BUILDTYPE - STATIC; SHARED; EXECUTABLE
## PUBLIC_LINK  - public libraries to link
## PRIVATE_LINK - private libraries to link
## FILES     - adds these files to the target
## FOLDER    - adds all files in these folders to the targets
## EXCLUDE   - excludes these files from the target
##
## note: If no files and folders are passed, all files from the level of current
##       CMakeLists.txt are added recursively.
##
#################################################################################
include(CMake/FileUtilities.cmake)

function(vf_add_library)

    set( options )
    set( oneValueArgs NAME BUILDTYPE)
    set( multiValueArgs PUBLIC_LINK PRIVATE_LINK FILES FOLDER EXCLUDE MODULEFOLDER)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    vf_add_target(NAME ${ARG_NAME} BUILDTYPE ${ARG_BUILDTYPE} PUBLIC_LINK ${ARG_PUBLIC_LINK} PRIVATE_LINK ${ARG_PRIVATE_LINK} FILES ${ARG_FILES} FOLDER ${ARG_FOLDER} EXCLUDE ${ARG_EXCLUDE} MODULEFOLDER ${ARG_MODULEFOLDER})

    if(NOT DEFINED ARG_BUILDTYPE)
        string(REGEX REPLACE "src" "tests/unit-tests" test_path "${CMAKE_CURRENT_SOURCE_DIR}")
        string(REGEX REPLACE "src" "tests/unit-tests" test_build_path "${CMAKE_CURRENT_BINARY_DIR}")
        if (EXISTS ${test_path})
            add_subdirectory(${test_path} ${test_build_path})
        endif()
    endif()
endfunction()

function(vf_add_tests2)

    if (NOT VF_ENABLE_UNIT_TESTS)
        return()
    endif()

    set( options )
    set( oneValueArgs NAME)
    set( multiValueArgs PUBLIC_LINK PRIVATE_LINK FILES FOLDER EXCLUDE MODULEFOLDER)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if(DEFINED ARG_NAME) 
        set(library_test_name "${ARG_NAME}")
    else()
        vf_get_library_test_name(library_test_name)
    endif()


    vf_add_target(NAME ${library_test_name} BUILDTYPE binary PUBLIC_LINK ${ARG_PUBLIC_LINK} PRIVATE_LINK ${ARG_PRIVATE_LINK} FILES ${ARG_FILES} FOLDER ${ARG_FOLDER} EXCLUDE ${ARG_EXCLUDE} MODULEFOLDER ${ARG_MODULEFOLDER})

    group_target (${library_test_name} ${testFolder})

    # link googlemock
    target_link_libraries(${library_test_name} PRIVATE GTest::gmock_main)

    # add the target to ctest
    gtest_add_tests(TARGET ${library_test_name})
    
endfunction()


function(vf_add_target)
    set( options )
    set( oneValueArgs NAME BUILDTYPE)
    set( multiValueArgs PUBLIC_LINK PRIVATE_LINK FILES FOLDER EXCLUDE MODULEFOLDER)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if(DEFINED ARG_NAME) 
        set(library_name ${ARG_NAME})
    else()
        vf_get_library_name (library_name)
    endif()
    vf_get_library_name (folder_name) # folder_name is not equal to library_name when ARG_NAME was set

    if(NOT DEFINED ARG_BUILDTYPE)
        if(BUILD_SHARED_LIBS)
            set(ARG_BUILDTYPE "shared")
        else()
            set(ARG_BUILDTYPE "static")
        endif()
    endif()

    if(DEFINED ARG_MODULEFOLDER)
        set(folder_name ${ARG_MODULEFOLDER})
    endif()

    # status("Configuring the target: ${library_name} (type=${ARG_BUILDTYPE})...")

    #################################################################
    ###   FIND FILES                                              ###
    #################################################################
    collectFiles(sourceFiles "${ARG_FILES}" "${ARG_FOLDER}" "${ARG_EXCLUDE}")
    includeAllFiles (${folder_name} "${sourceFiles}")

    #################################################################
    ###   ADD TARGET                                              ###
    #################################################################
    if(${ARG_BUILDTYPE} MATCHES binary)
        add_executable(${library_name} ${MY_SRCS} )
        group_target(${library_name} ${appFolder})
    elseif(${ARG_BUILDTYPE} MATCHES shared)
        add_library(${library_name} SHARED ${MY_SRCS} )
        group_target(${library_name} ${libraryFolder})
        elseif(${ARG_BUILDTYPE} MATCHES static)
        add_library(${library_name} STATIC ${MY_SRCS} )
        group_target(${library_name} ${libraryFolder})
    else()
        message(FATAL_ERROR "build_type=${ARG_BUILDTYPE} doesn't match BINARY, SHARED or STATIC")
    endif()

    set_target_properties(${library_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

    #################################################################
    ###   ADDITIONAL LINK LIBRARIES                               ###
    #################################################################
    status_lib("Link Depending public libraries: ${ARG_PUBLIC_LINK}")
    status_lib("Link Depending private libraries: ${ARG_PRIVATE_LINK}")
    if (ARG_PUBLIC_LINK)
        target_link_libraries(${library_name} PUBLIC ${ARG_PUBLIC_LINK})
    endif()
    if (ARG_PRIVATE_LINK)
        target_link_libraries(${library_name} PRIVATE ${ARG_PRIVATE_LINK})
    endif()

    status("Target: ${library_name} (type=${ARG_BUILDTYPE}) configured.")


endfunction()

#################################################################################
## Add a test executable corresponding to the added target.
## Must be called after vf_add_library().
## The name of the test executable is: vf_get_library_name()Tests
##
## Precondition: VF_ENABLE_UNIT_TESTS needs to be ON
#################################################################################
function(vf_add_tests)

    if (NOT VF_ENABLE_UNIT_TESTS)
        return()
    endif()

    set( options )
    set( oneValueArgs NAME)
    set( multiValueArgs)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if(DEFINED ARG_NAME) 
        set(library_test_name "${ARG_NAME}Tests")
        set(library_name "${ARG_NAME}")
    else()
        vf_get_library_test_name(library_test_name)
        vf_get_library_name (library_name)
    endif()

    vf_get_library_name (folder_name)

    status("Configuring test executable: ${library_test_name}")

    # set test files to MY_SRCS
    file ( GLOB_RECURSE all_files ${VIRTUAL_FLUIDS_GLOB_FILES} )
    includeTestFiles (${folder_name} "${all_files}")

    # add the target
    add_executable(${library_test_name} ${MY_SRCS})
    group_target (${library_test_name} ${testFolder})

    set_target_properties(${library_test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


    # link tested library
    target_link_libraries(${library_test_name} PRIVATE ${library_name} project_options)

    # link tested library
    target_include_directories(${library_test_name} PRIVATE ${CMAKE_BINARY_DIR})
    target_include_directories(${library_test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${library_test_name} PRIVATE ${VF_SRC_DIR})

    # link googlemock
    target_link_libraries(${library_test_name} PRIVATE GTest::gmock_main)

    # add the target to ctest
    gtest_add_tests(TARGET ${library_test_name})

endfunction()

#################################################################################
## group target for VisualStudio
#################################################################################
function(group_target targetName folderName)
    set_property( TARGET  ${targetName}  PROPERTY  FOLDER  ${folderName} )
endfunction(group_target)

#################################################################################
## load user apps, which are specified in the machine file
#################################################################################
function(vf_load_user_apps)
    foreach(app IN LISTS USER_APPS)
      add_subdirectory(${app})
    endforeach()
endfunction()
