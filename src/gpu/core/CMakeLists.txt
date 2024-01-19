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
#  You should have received a copy of the GNU General Public License along
#  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
#
#  SPDX-License-Identifier: GPL-3.0-or-later
#  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
#
#! \author Soeren Peters
#=======================================================================================
project(gpu_core LANGUAGES CUDA CXX)


vf_add_library(NAME gpu_core PUBLIC_LINK basics lbm parallel logger project_options project_warnings PRIVATE_LINK GridGenerator cuda_helper)

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

set_target_properties(gpu_core PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)

target_include_directories(gpu_core SYSTEM PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} "${VF_THIRD_DIR}/cuda_samples/")

if(MSVC)
    set_target_properties(gpu_core PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

## WARNINGS
# TODO: Remove this and fix those warnings: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/163
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    list(APPEND WARNINGS "-Wold-style-cast;-Wno-old-style-cast")
    list(APPEND WARNINGS "-Wsign-conversion;-Wno-sign-conversion")
    list(APPEND WARNINGS "-Wdouble-promotion;-Wno-double-promotion")
    list(APPEND WARNINGS "-Wshadow;-Wno-shadow")

    list(APPEND WARNINGS "-Wextra") # Wno-unused-parameter comes with -Wextra
    list(APPEND WARNINGS "-Wall") # Wno-unused-function comes with -Wall
    list(APPEND WARNINGS "-Wunused") # -Wno-unused-function comes with -Wunused
    list(APPEND WARNINGS "-Wno-unused-parameter")
    list(APPEND WARNINGS "-Wno-unused-function")
    list(APPEND WARNINGS "-Wno-implicit-float-conversion")
    list(APPEND WARNINGS "-Wno-comment")
    list(APPEND WARNINGS "-Wno-reorder-ctor")
    list(APPEND WARNINGS "-Wno-uninitialized")
    list(APPEND WARNINGS "-Wno-reorder")
    list(APPEND WARNINGS "-Wno-strict-aliasing")

    if(VF_ENABLE_UNIT_TESTS)
        target_compile_options(gpu_coreTests PUBLIC "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wnull-dereference;-Xcompiler=-Wno-null-dereference>")
        target_compile_options(gpu_coreTests PUBLIC "$<$<COMPILE_LANGUAGE:CXX>:-Wnull-dereference;-Wno-null-dereference>")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    list(APPEND WARNINGS "-Wuseless-cast;-Wno-useless-cast")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    list(APPEND WARNINGS "/wd4244") # conversion
    list(APPEND WARNINGS "/wd4100") # unreferenced formal parameter
    list(APPEND WARNINGS "/wd4701") # potentially uninitialized local variable
    list(APPEND WARNINGS "/wd4703") # potentially uninitialized local pointer
    list(APPEND WARNINGS "/wd4305") # 'initializing': truncation from double to float  
    list(APPEND WARNINGS "/wd4505") # 'XXX': unreferenced local function has been removed
    list(APPEND WARNINGS "/wd4324") # structure was padded due to alignment specifier
    list(APPEND WARNINGS "/w14265;/wd4265")
    list(APPEND WARNINGS "/wd4389")
    list(APPEND WARNINGS "/w14296;/wd4296")

    target_link_options(gpu_core PUBLIC $<DEVICE_LINK:-Xcompiler=/wd4100>)
endif()

foreach(flag IN LISTS WARNINGS)
    target_compile_options(gpu_core PUBLIC "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${flag}>")
    target_compile_options(gpu_core PUBLIC "$<$<COMPILE_LANGUAGE:CXX>:${flag}>")
endforeach()
