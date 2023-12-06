# #######################################################################################
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
#  You should have received a copy of the GNU General Public License along
#  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
# #######################################################################################
include(CMake/compilerflags/CompilerWarnings.cmake)
include(CMake/compilerflags/CompilerOptions.cmake)
include(CMakePrintHelpers)
include(CMake/VirtualFluidsMacros.cmake)

# prevent in-source builds
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds are prohibited. Create a new directory and build there.")
endif ()

# check if we build static
if(${BUILD_SHARED_LIBS})
    message(FATAL_ERROR "VirtualFluids only supports static build. Rerun cmake with -DBUILD_SHARED_LIBS=OFF")
endif()

# set build type to Release when not set
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# global cmake options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ".cmake")
set(libraryFolder "libs")
set(testFolder    "tests")
set(appFolder     "apps")
set(thirdFolder   "3rd")

set(VF_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set(VF_THIRD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty)
set(VF_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(VF_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# windows: use multi-threaded dynamically-linked runtime library
if(BUILD_SHARED_LIBS)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

# bindings
if (VF_ENABLE_PYTHON_BINDINGS)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# project warnings and project options
add_library(project_warnings INTERFACE)
add_library(project_options INTERFACE)


set_project_warnings(project_warnings)
set_project_options(project_options)

include(CMake/Sanitizers.cmake)
enable_sanitizers(project_options)

if(VF_ENABLE_CACHE)
    include(CMake/Cache.cmake)
    enable_cache()
endif()

include(CMake/StaticAnalyzers.cmake)

if(VF_ENABLE_CLANG_TIDY)
    enable_clang_tidy(project_options)
endif()

if(VF_ENABLE_CPPCHECK)
    enable_cppcheck()
endif()

if(VF_ENABLE_INCLUDE_WHAT_YOU_USE)
    enable_include_what_you_use()
endif()

# set gpu features
if(VF_ENABLE_GPU)
    include(CheckLanguage)
    check_language(CUDA)

    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR "CUDA Compiler was requested but is not found on the system.")
    endif()

    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)

    enable_language(CUDA)

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        message(WARNING "CMAKE_CUDA_ARCHITECTURES was not defined and is set to 30 (CUDA support until 10.1 only).")
        set(CMAKE_CUDA_ARCHITECTURES 30)
    endif()

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")

    message(STATUS "CMAKE_CUDA_FLAGS: ${CMAKE_CUDA_FLAGS}")
    message(STATUS "CUDA Architecture: ${CMAKE_CUDA_ARCHITECTURES}")
    set(CMAKE_CUDA_ARCHITECTURES
        "${CMAKE_CUDA_ARCHITECTURES}"
        CACHE STRING "Cuda Architecture (compute capabilitiy)")

    set(CMAKE_CUDA_FLAGS_DEBUG
        " -G"
        CACHE STRING "" FORCE)

    # we disable the usage of cuda response files here usually CUDA_INCLUDES.rsp is genereated by cmake containing all
    # include paths and is passed in compile_commands.json via the --options-file flag this .rsp file can not be parsed
    # by clangd and therefore we disable it
    set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_INCLUDES 0)
endif()


#################################################################
### load machine file                                         ###
#################################################################
site_name(MACHINE_NAME)
string(TOUPPER  "${MACHINE_NAME}" MACHINE_NAME)

set(MACHINE_FILE "${CMAKE_CURRENT_LIST_DIR}/cmake_config_files/${MACHINE_NAME}.config.cmake")

IF(NOT EXISTS ${MACHINE_FILE})
    status("No configuration file found: ${MACHINE_FILE}.")
ELSE()
    status("Load configuration file: ${MACHINE_FILE}")
    include(${MACHINE_FILE})
ENDIF()


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
