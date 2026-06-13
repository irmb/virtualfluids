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
#  SPDX-License-Identifier: GPL-3.0-or-later
#  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
#
# #################################################################################

if(VF_ENABLE_GPU)
    set(PYFLUIDS_LANGUAGES CUDA CXX)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF) # need it to work with cuda
else()
    set(PYFLUIDS_LANGUAGES CXX)
endif()

project(VirtualFluidsPython LANGUAGES ${PYFLUIDS_LANGUAGES})

set(PYFLUIDS_DIR ${SKBUILD_PLATLIB_DIR}/pyfluids)
set(PYMUPARSER_DIR ${SKBUILD_PLATLIB_DIR}/pymuparser)

pybind11_add_module(python_bindings pyfluids.cpp)
target_link_libraries(python_bindings PRIVATE logger basics parallel lbm project_options)

set_target_properties(  python_bindings PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                        OUTPUT_NAME "_pyfluids"
                        POSITION_INDEPENDENT_CODE ON)

IF(VF_ENABLE_GPU)
    target_compile_definitions(python_bindings PUBLIC VF_ENABLE_GPU)
    target_include_directories(python_bindings PRIVATE ${VF_THIRD_DIR}/cuda_samples/ )
    target_link_libraries(python_bindings PRIVATE GridGenerator gpu_core)
ENDIF()

IF(VF_ENABLE_CPU)
    target_compile_definitions(python_bindings PUBLIC VF_ENABLE_CPU)
    target_link_libraries(python_bindings PRIVATE simulationconfig cpu_core muparser)

    # copied from src/cpu/core/CMakeLists.txt
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
        # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
        target_compile_options(python_bindings PUBLIC "-Wextra") # Wno-unused-parameter comes with -Wextra
        target_compile_options(python_bindings PUBLIC "-Wall") # Wno-unused-function comes with -Wall
        target_compile_options(python_bindings PUBLIC "-Wunused") # -Wno-unused-function comes with -Wunused
        target_compile_options(python_bindings PUBLIC "-Wno-unused-parameter")
        target_compile_options(python_bindings PUBLIC "-Wno-unused-function")
        target_compile_options(python_bindings PUBLIC "-Wno-reorder")
        
        # deduplicate option to silence warning about null dereference in interactor3d
        target_compile_options(python_bindings PRIVATE "-Wnull-dereference")
        target_compile_options(python_bindings PRIVATE "-Wno-null-dereference")
    endif()
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        target_compile_options(${library_name} PUBLIC "/wd4244") # conversion
        target_compile_options(${library_name} PUBLIC "/wd4100") # unreferenced formal parameter
        target_compile_options(${library_name} PUBLIC "/wd4701") # potentially uninitialized local variable
        target_compile_options(${library_name} PUBLIC "/wd4703") # potentially uninitialized local pointer
        target_compile_options(${library_name} PUBLIC "/wd4305") #'initializing': truncation from double to float
        target_compile_options(${library_name} PUBLIC "/wd4505") # C4505: 'XXX': unreferenced local function has been removed
        target_compile_options(${library_name} PUBLIC "/wd4310") # C4310: cast truncates constant value, triggerd by muParserbase.h
        target_compile_options(${library_name} PUBLIC "/wd4456") # C4456: declaration of hides previous local declaration
    endif()
    # include bindings for muparsers
    pybind11_add_module(pymuparser MODULE muParser.cpp)

    # TODO: Move this to MuParser CMakeLists.txt
    set_target_properties(muparser PROPERTIES POSITION_INDEPENDENT_CODE ON)

    set_target_properties(  pymuparser PROPERTIES
                            LIBRARY_OUTPUT_DIRECTORY ${PYMUPARSER_DIR}
                            OUTPUT_NAME "bindings")
    target_link_libraries(pymuparser PRIVATE muparser project_options)
    add_dependencies(python_bindings pymuparser)
ENDIF()
