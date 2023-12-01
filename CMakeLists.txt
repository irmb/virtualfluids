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

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(
    VirtualFluids
    VERSION 0.1.0
    DESCRIPTION "CFD code based on the Lattice Boltzmann Method"
    HOMEPAGE_URL "https://www.tu-braunschweig.de/irmb/forschung/virtualfluids"
    LANGUAGES CXX)

#############################################################
###                    Options                            ###
#############################################################
option(BUILD_SHARED_LIBS "Building VirtualFluids with static libs. Shared libs are not supported." OFF)

option(VF_ENABLE_CPU "Build VirtualFluids cpu variant" OFF)
option(VF_ENABLE_GPU "Build VirtualFluids gpu variant" OFF)
option(VF_ENABLE_ALL_APPS "Build All example applications" OFF)
option(VF_ENABLE_UNIT_TESTS "Build VirtualFluids unit tests" OFF)
option(VF_ENABLE_PYTHON_BINDINGS "Build the python bindings" OFF)
option(VF_ENABLE_DOUBLE_ACCURACY "Use double accuracy" OFF)

# --- TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/139 ---
# option(VF_ENABLE_MPI "include MPI library support" ON)
option(VF_ENABLE_OPENMP "Build VirtualFluids with openmp" ON)
option(VF_ENABLE_BOOST "Build VirtualFluids with boost" OFF)

# analyzer
option(VF_ENABLE_WARNINGS_AS_ERRORS "Make warnings to errors to identify problems fast" OFF)
option(VF_ENABLE_CLANG_TIDY "Add the clang tidy checks to the targets" OFF)
option(VF_ENABLE_INCLUDE_WHAT_YOU_USE "Add IWYU to the targets" OFF)
option(VF_ENABLE_CPPCHECK "Add cppcheck to the targets" OFF)
option(VF_ENABLE_CACHE "Enable ccache" ON)

option(VF_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
option(VF_ENABLE_SANITIZER_UNDEFINED_BEHAVIOR "Enable undefined behavior sanitizer" OFF)
option(VF_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
option(VF_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
option(VF_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)

#############################################################
###                  Libraries                            ###
#############################################################
include(CMake/ProjectSettings.cmake)
include(CMake/Packages.cmake)


add_subdirectory(src/logger)
add_subdirectory(src/basics)
add_subdirectory(src/parallel)
add_subdirectory(src/lbm)


if(VF_ENABLE_CPU)
    include(cpu.cmake)
endif()

if(VF_ENABLE_GPU)
    include(gpu.cmake)
endif()

if(VF_ENABLE_PYTHON_BINDINGS)
    add_subdirectory(pythonbindings)
endif()

vf_load_user_apps()
