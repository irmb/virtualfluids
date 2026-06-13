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
#  SPDX-FileCopyrightText: Copyright �� VirtualFluids Project contributors, see AUTHORS.md in root folder
#
# \author Hussein Alihussein
#  This file was partially implemented with the help of an AI model (OpenAI ChatGPT; specific model version unknown).
#=======================================================================================
project(ResearchCrossingNeighbourhoodGpu)

vf_add_executable(NAME ResearchCrossingNeighbourhoodGpu PUBLIC_LINK gpu_core PRIVATE_LINK gpu_core)

# Fast Winding (libigl) for the GPU app
set(VF_FAST_WINDING_ROOT "" CACHE PATH "Path to fast-winding-number-soups checkout")
if(NOT VF_FAST_WINDING_ROOT)
  set(_vf_fast_default "${CMAKE_SOURCE_DIR}/3rdParty/fast-winding-number-soups")
  if(EXISTS "${_vf_fast_default}/libigl/include/igl/FastWindingNumberForSoups.h")
    set(VF_FAST_WINDING_ROOT "${_vf_fast_default}" CACHE PATH "Path to fast-winding-number-soups checkout" FORCE)
  endif()
endif()

if(VF_FAST_WINDING_ROOT)
  target_include_directories(ResearchCrossingNeighbourhoodGpu PRIVATE
      ${VF_FAST_WINDING_ROOT}/include
      ${VF_FAST_WINDING_ROOT}/libigl/include
  )
  target_compile_definitions(ResearchCrossingNeighbourhoodGpu PRIVATE VF_HAS_FAST_WINDING=1)
endif()

target_include_directories(ResearchCrossingNeighbourhoodGpu PRIVATE
    ${CMAKE_SOURCE_DIR}/src/cpu/core
    ${CMAKE_SOURCE_DIR}/src/cpu/core/LBM
    ${CMAKE_SOURCE_DIR}/src/cpu/core/BoundaryConditions
    ${CMAKE_SOURCE_DIR}/src/cpu/core/Simulation
    ${CMAKE_SOURCE_DIR}/src/cpu/core/Visitors
    ${CMAKE_SOURCE_DIR}/src/cpu/core/SimulationObservers
    ${CMAKE_SOURCE_DIR}/src/cpu/core/Interactors
    ${CMAKE_SOURCE_DIR}/src/cpu/SIMULATIONCONFIG
    ${CMAKE_SOURCE_DIR}/3rdParty/MuParser/include)

target_link_directories(ResearchCrossingNeighbourhoodGpu PRIVATE
    ${CMAKE_BINARY_DIR}/lib)

# Silence third-party libigl fast-winding warnings promoted to errors
if(NOT MSVC)
  target_compile_options(ResearchCrossingNeighbourhoodGpu PRIVATE
      -Wno-error=pedantic
      -Wno-error=sign-compare
      -Wno-error=type-limits
      -Wno-error=deprecated-declarations
      -Wno-error=unused-variable
      -Wno-error=cast-align
      $<$<CXX_COMPILER_ID:GNU>:-Wno-error=class-memaccess>)
else()
  target_compile_options(ResearchCrossingNeighbourhoodGpu PRIVATE
      /wd4201 # nameless struct/union in libigl fast-winding
      /wd4267 # size_t to int conversions in libigl fast-winding
      /wd4456 # unused lambda capture in libigl fast-winding
      /wd4189 # unused locals in libigl fast-winding
  )
endif()

if(VF_ENABLE_CPU)
    target_link_libraries(ResearchCrossingNeighbourhoodGpu PRIVATE
        cpu_core
        muparser)
endif()
