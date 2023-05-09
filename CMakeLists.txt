#################################################################################
#   _    ___      __              __________      _     __
# | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
# | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
# | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
# |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
#
#################################################################################
cmake_minimum_required(VERSION 3.15..3.20 FATAL_ERROR)

project(VirtualFluids
        VERSION 0.1.0
        DESCRIPTION "CFD code based on the Lattice Boltzmann Method"
        HOMEPAGE_URL "https://www.tu-braunschweig.de/irmb/forschung/virtualfluids"
        LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ".cmake")
set(libraryFolder "libs")
set(testFolder    "tests")
set(appFolder     "apps")
set(thirdFolder   "3rd")

set (VF_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set (VF_THIRD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty)
set (VF_SRC_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/src)
set (VF_ROOT_DIR  ${CMAKE_CURRENT_SOURCE_DIR})

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are prohibited. "
    "Create a new directory and build there.")
endif ()

#################################################################################
#  OPTIONS
#################################################################################
option(BUILD_VF_CPU "Build VirtualFluids cpu variant" OFF)
option(BUILD_VF_GPU "Build VirtualFluids gpu variant" OFF)

option(BUILD_USE_OPENMP "Build VirtualFluids with openmp" ON)
option(BUILD_USE_BOOST "Build VirtualFluids with boost" OFF)
option(BUILD_USE_MPI "include MPI library support" ON)

# vf gpu
option(BUILD_NUMERIC_TESTS   "Build numeric tests"         OFF)

option(BUILD_VF_UNIT_TESTS "Build VirtualFluids unit tests" OFF)
option(BUILD_VF_CLANG_TIDY "Add the clang tidy checks to the targets" OFF)
option(BUILD_VF_INCLUDE_WHAT_YOU_USE "Add IWYU to the targets" OFF)
option(BUILD_VF_CPPCHECK "Add cppcheck to the targets" OFF)
option(BUILD_VF_COVERAGE "Add the -coverage compiler flag." OFF)

option(BUILD_CUDA_LTO "Enables the cuda link optimization." OFF)

option(BUILD_SHARED_LIBS "" OFF)
option(BUILD_WARNINGS_AS_ERRORS "" OFF)

# windows: use multi-threaded dynamically-linked runtime library
if(BUILD_SHARED_LIBS)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

option(BUILD_VF_PYTHON_BINDINGS "" OFF)

option(BUILD_VF_DOUBLE_ACCURACY "Use double accuracy" OFF)
option(BUILD_VF_ALL_SAMPLES "Build All Samples" OFF)

#################################################################################
#  MACROS
#################################################################################
add_library(project_warnings INTERFACE)

# standard compiler warnings - Link this 'library' to use the warnings specified in CompilerWarnings.cmake
include(${VF_CMAKE_DIR}/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# enable cache system
include(${VF_CMAKE_DIR}/Cache.cmake)

include(CMakePrintHelpers)
include(${VF_CMAKE_DIR}/VirtualFluidsMacros.cmake)

IF( BUILD_VF_DOUBLE_ACCURACY )
    list(APPEND VF_COMPILER_DEFINITION VF_DOUBLE_ACCURACY)
ENDIF()

# set gpu features
if(BUILD_VF_GPU)
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


    if(BUILD_CUDA_LTO)
       if(CMAKE_CUDA_ARCHITECTURES LESS 50)
            message(FATAL_ERROR "CUDA Link time optimization requires CUDA 11.2 and CC 5.0.")
       else()

        set(CMAKE_CUDA_FLAGS "-dlto -arch=sm_${CMAKE_CUDA_ARCHITECTURES}")
        set(CMAKE_CUDA_ARCHITECTURES OFF)

        list(APPEND VF_COMPILER_DEFINITION BUILD_CUDA_LTO)

       endif()
    endif()

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")

    message(STATUS "CMAKE_CUDA_FLAGS: ${CMAKE_CUDA_FLAGS}")
    message(STATUS "CUDA Architecture: ${CMAKE_CUDA_ARCHITECTURES}")
    set(CMAKE_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}" CACHE STRING "Cuda Architecture (compute capabilitiy)")

    set(CMAKE_CUDA_FLAGS_DEBUG " -G" CACHE STRING "" FORCE)
endif()


#################################################################################
#  COMMON LIBRARIES
#################################################################################
include(FetchContent)

# spdlog
set(spdlog_version "v1.9.1")
set(spdlog_url "https://github.com/gabime/spdlog")
message(STATUS "Fetching spdlog: ${spdlog_version}")
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY ${spdlog_url}
        GIT_TAG        ${spdlog_version}
)

FetchContent_MakeAvailable(spdlog)
target_compile_options(spdlog PRIVATE "-fPIC")
groupTarget(spdlog ${thirdFolder})

# googletest
if(BUILD_VF_UNIT_TESTS)
    FetchContent_Declare(
        googletest
        DOWNLOAD_EXTRACT_TIMESTAMP FALSE # https://cmake.org/cmake/help/latest/policy/CMP0135.html
        URL https://github.com/google/googletest/archive/1f643f71d4151c3b364c0e9302042f7a6debd439.zip # 30.11.2022
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(googletest)

    groupTarget(gmock ${thirdFolder}/googletest)
    groupTarget(gmock_main ${thirdFolder}/googletest)
    groupTarget(gtest ${thirdFolder}/googletest)
    groupTarget(gtest_main ${thirdFolder}/googletest)

    include(GoogleTest)
    enable_testing()
endif()

if(BUILD_USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

if(BUILD_USE_MPI)
    find_package(MPI REQUIRED)
    list(APPEND VF_COMPILER_DEFINITION VF_MPI)
endif()

# boost
IF(BUILD_USE_BOOST)
    list(APPEND VF_COMPILER_DEFINITION BUILD_USE_BOOST)

    set(Boost_USE_STATIC_LIBS ON)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_STATIC_RUNTIME ON)

    # minimum boost version: 1.60
    # no packages specfied - only headeronly libraries
    find_package(Boost 1.60 REQUIRED)
ENDIF()

add_subdirectory(src/logger)
add_subdirectory(src/basics)
add_subdirectory(src/mpi)
add_subdirectory(src/lbm)


#################################################################################
#  VIRTUAL FLUIDS CPU / GPU
#################################################################################
if (BUILD_VF_CPU)
    include(cpu.cmake)
endif()
if(BUILD_VF_GPU)
    add_subdirectory(src/cuda)
    include(gpu.cmake)
endif()

if (BUILD_VF_PYTHON_BINDINGS)
    set(pybind_version "v2.10.4")
    set(pybind_url "https://github.com/pybind/pybind11")
    message(STATUS "Fetching pybind: ${pybind_version}")
    FetchContent_Declare(
            pybind11
            GIT_REPOSITORY ${pybind_url}
            GIT_TAG        ${pybind_version}
            OVERRIDE_FIND_PACKAGE
    )

    find_package(pybind11)
    add_subdirectory(pythonbindings)
endif()


vf_load_user_apps()
