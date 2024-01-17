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
vf_add_library(NAME GridGenerator PUBLIC_LINK basics logger OpenMP::OpenMP_CXX PRIVATE_LINK project_options project_warnings)

target_include_directories(GridGenerator PRIVATE "${VF_ROOT_DIR}/src/gpu/GridGenerator/")


## WARNINGS
# TODO: Remove this and fix those warnings: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/163
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    target_compile_options(GridGenerator PUBLIC "-Wold-style-cast;-Wno-old-style-cast")
    target_compile_options(GridGenerator PUBLIC "-Wsign-conversion;-Wno-sign-conversion")
    target_compile_options(GridGenerator PUBLIC "-Wdouble-promotion;-Wno-double-promotion")
    target_compile_options(GridGenerator PUBLIC "-Wshadow;-Wno-shadow")
    target_compile_options(GridGenerator PRIVATE "-Wcast-align;-Wno-cast-align")

    target_compile_options(GridGenerator PRIVATE "-Wextra") # Wno-unused-parameter comes with -Wextra
    target_compile_options(GridGenerator PRIVATE "-Wall") # Wno-unused-function comes with -Wall
    target_compile_options(GridGenerator PRIVATE "-Wunused") # -Wno-unused-function comes with -Wunused
    target_compile_options(GridGenerator PRIVATE "-Wno-unused-parameter")
    target_compile_options(GridGenerator PRIVATE "-Wno-unused-function")
    target_compile_options(GridGenerator PRIVATE "-Wno-implicit-float-conversion")
    target_compile_options(GridGenerator PRIVATE "-Wno-comment")
    target_compile_options(GridGenerator PRIVATE "-Wno-reorder-ctor")
    target_compile_options(GridGenerator PRIVATE "-Wno-uninitialized")
    target_compile_options(GridGenerator PRIVATE "-Wno-reorder")
    target_compile_options(GridGenerator PRIVATE "-Wno-strict-aliasing")
    target_compile_options(GridGenerator PRIVATE "-Wformat=2;-Wno-format-nonliteral")
    target_compile_options(GridGenerator PRIVATE "-Wno-null-dereference")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(GridGenerator PUBLIC "-Wuseless-cast;-Wno-useless-cast")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(GridGenerator PRIVATE "/wd4244") # conversion
    target_compile_options(GridGenerator PRIVATE "/wd4100") # unreferenced formal parameter
    target_compile_options(GridGenerator PRIVATE "/wd4701") # potentially uninitialized local variable
    target_compile_options(GridGenerator PRIVATE "/wd4703") # potentially uninitialized local pointer
    target_compile_options(GridGenerator PRIVATE "/wd4305") #'initializing': truncation from double to float
    target_compile_options(GridGenerator PRIVATE "/wd4389") #''==': signed/unsigned mismatch
endif()