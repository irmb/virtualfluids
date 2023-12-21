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

set(VF_LIBRARIES)

IF(METIS_RELEASE_LIBRARY AND METIS_DEBUG_LIBRARY)
   list(APPEND VF_LIBRARIES optimized ${METIS_RELEASE_LIBRARY} debug ${METIS_DEBUG_LIBRARY})
ELSE()
   list(APPEND VF_LIBRARIES metis)
ENDIF()

IF(${VF_CPU_ENABLE_VTK})
   list(APPEND VF_LIBRARIES optimized ${VTK_LIBRARIES} debug ${VTK_LIBRARIES})
ENDIF()

IF(${VF_CP_ENABLE_CATALYST})
   list(APPEND VF_LIBRARIES optimized vtkPVPythonCatalyst debug vtkPVPythonCatalyst)
   list(APPEND VF_LIBRARIES optimized vtkParallelMPI debug vtkParallelMPI )
ENDIF()

if(VF_ENABLE_OPENMP)
    list(APPEND VF_LIBRARIES OpenMP::OpenMP_CXX)
endif()

set(library_name cpu_core)

vf_add_library(NAME ${library_name} PUBLIC_LINK basics muparser ${VF_LIBRARIES} parallel project_options PRIVATE_LINK lbm logger project_warnings)

target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Connectors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Data)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Interactors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LBM)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LBM/Interpolation)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Simulation)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Visitors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/SimulationObservers)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Utilities)


IF(${VF_CPU_ENABLE_VTK})
   target_include_directories(${library_name} PRIVATE ${VTK_INCLUDE_DIRS})
ENDIF()

# TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/139
# if(VF_ENABLE_MPI)
    target_link_libraries(${library_name} PRIVATE MPI::MPI_CXX)
# endif()

## WARNINGS
# TODO: Remove this and fix those warnings: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/163
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    target_compile_options(${library_name} PUBLIC "-Wextra") # Wno-unused-parameter comes with -Wextra
    target_compile_options(${library_name} PUBLIC "-Wall") # Wno-unused-function comes with -Wall
    target_compile_options(${library_name} PUBLIC "-Wunused") # -Wno-unused-function comes with -Wunused
    target_compile_options(${library_name} PUBLIC "-Wno-unused-parameter")
    target_compile_options(${library_name} PUBLIC "-Wno-unused-function")
    target_compile_options(${library_name} PUBLIC "-Wno-reorder")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(${library_name} PUBLIC "/wd4244") # conversion
    target_compile_options(${library_name} PUBLIC "/wd4100") # unreferenced formal parameter
    target_compile_options(${library_name} PUBLIC "/wd4701") # potentially uninitialized local variable
    target_compile_options(${library_name} PUBLIC "/wd4703") # potentially uninitialized local pointer
    target_compile_options(${library_name} PUBLIC "/wd4458") # declaration of 'XXX' hides class member
    target_compile_options(${library_name} PUBLIC "/wd4305") #'initializing': truncation from double to float
    target_compile_options(${library_name} PUBLIC "/wd4505") # C4505: 'XXX': unreferenced local function has been removed
    target_compile_options(${library_name} PUBLIC "/wd4310") # C4310: cast truncates constant value, triggerd by muParserbase.h
    target_compile_options(${library_name} PUBLIC "/wd4005") # disable macro redefinition (triggered by metis.h)
    target_compile_options(${library_name} PUBLIC "/wd4456") # C4456: declaration of hides previous local declaration
endif()