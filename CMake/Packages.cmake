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
include(FetchContent)

# spdlog
set(spdlog_version "v1.14.1")
set(spdlog_url "https://github.com/gabime/spdlog")
message(STATUS "Fetching spdlog: ${spdlog_version}")
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY ${spdlog_url}
    GIT_TAG ${spdlog_version}
    GIT_SHALLOW 1)

FetchContent_MakeAvailable(spdlog)
if(NOT MSVC)
    target_compile_options(spdlog PRIVATE "-fPIC")
endif()
if(MSVC)
    target_compile_options(spdlog PUBLIC "$<$<COMPILE_LANGUAGE:CXX>:/wd4996>")
    target_compile_options(spdlog PUBLIC "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/wd4996>")
endif()
group_target(spdlog ${thirdFolder})

# googletest
if(VF_ENABLE_UNIT_TESTS)
    FetchContent_Declare(
        googletest
        DOWNLOAD_EXTRACT_TIMESTAMP FALSE # https://cmake.org/cmake/help/latest/policy/CMP0135.html
        URL https://github.com/google/googletest/archive/1f643f71d4151c3b364c0e9302042f7a6debd439.zip # 30.11.2022
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt
        ON
        CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(googletest)

    group_target(gmock ${thirdFolder}/googletest)
    group_target(gmock_main ${thirdFolder}/googletest)
    group_target(gtest ${thirdFolder}/googletest)
    group_target(gtest_main ${thirdFolder}/googletest)

    if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
        target_compile_options(gtest PRIVATE "-Wno-implicit-int-float-conversion")
    endif()

    include(GoogleTest)
    enable_testing()
endif()

if(VF_ENABLE_OPENMP)
    # Only require OpenMP for C++ to avoid attempting OpenMP offloading for CUDA,
    # which is typically unsupported with NVCC on Windows.
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    target_compile_definitions(project_options INTERFACE VF_OPENMP)
endif()

# TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/139 if(VF_ENABLE_MPI)
find_package(MPI REQUIRED)
target_compile_definitions(project_options INTERFACE VF_MPI)
# endif()

# boost
if(VF_ENABLE_BOOST)
    target_compile_definitions(project_options INTERFACE VF_BOOST)

    set(Boost_USE_STATIC_LIBS ON)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_STATIC_RUNTIME ON)

    # minimum boost version: 1.60 no packages specfied - only headeronly libraries
    find_package(Boost 1.60 REQUIRED)
endif()

if(VF_ENABLE_PYTHON_BINDINGS)
    set(pybind_version "v3.0.2")
    set(pybind_url "https://github.com/pybind/pybind11")
    message(STATUS "Fetching pybind: ${pybind_version}")
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY ${pybind_url}
        GIT_TAG ${pybind_version}
        GIT_SHALLOW 1)

    FetchContent_MakeAvailable(pybind11)
endif()

# YAML
set(yaml_git_version "yaml-cpp-0.9.0")
set(yaml_url "https://github.com/jbeder/yaml-cpp")
message(STATUS "Fetching yaml-cpp: ${yaml_git_version}")
FetchContent_Declare(yaml-cpp
  GIT_REPOSITORY ${yaml_url}
  GIT_TAG ${yaml_git_version}
  GIT_SHALLOW 1) 
FetchContent_GetProperties(yaml-cpp)
if(NOT yaml-cpp_POPULATED)
  FetchContent_MakeAvailable(yaml-cpp)
endif()

# Metis
add_subdirectory(${VF_THIRD_DIR}/metis/metis-5.1.0)
target_compile_definitions(project_options INTERFACE VF_METIS)

# Fast Winding (libigl)
set(VF_FAST_WINDING_ROOT "" CACHE PATH "Path to fast-winding-number-soups checkout")
if(NOT VF_FAST_WINDING_ROOT)
    set(_vf_fast_default "${CMAKE_SOURCE_DIR}/3rdParty/fast-winding-number-soups")
    if(EXISTS "${_vf_fast_default}/libigl/include/igl/FastWindingNumberForSoups.h")
        set(VF_FAST_WINDING_ROOT "${_vf_fast_default}" CACHE PATH "Path to fast-winding-number-soups checkout" FORCE)
    endif()
endif()

function(vf_enable_fast_winding target_name)
    set(oneValueArgs VISIBILITY)
    cmake_parse_arguments(ARG "" "${oneValueArgs}" "" ${ARGN})

    if(NOT ARG_VISIBILITY)
        set(ARG_VISIBILITY PUBLIC)
    endif()

    if(NOT VF_FAST_WINDING_ROOT)
        return()
    endif()

    if(NOT EXISTS "${VF_FAST_WINDING_ROOT}/libigl/include/igl/FastWindingNumberForSoups.h")
        message(WARNING "VF_FAST_WINDING_ROOT does not contain libigl FastWindingNumberForSoups.h: ${VF_FAST_WINDING_ROOT}")
        return()
    endif()

    target_include_directories(${target_name} SYSTEM ${ARG_VISIBILITY}
        ${VF_FAST_WINDING_ROOT}/include
        ${VF_FAST_WINDING_ROOT}/libigl/include
    )
    target_compile_definitions(${target_name} ${ARG_VISIBILITY} VF_HAS_FAST_WINDING=1)
endfunction()
