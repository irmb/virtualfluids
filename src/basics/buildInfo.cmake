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
set(buildInfoPath ${CMAKE_BINARY_DIR}/buildInfo)
set(buildInfoFile buildInfo.cpp)
set(buildInfoInput ${CMAKE_CURRENT_LIST_DIR}/buildInfo.in.cpp)

include(${VF_CMAKE_DIR}/3rd/git/GetGitRevisionDescription.cmake)
get_git_head_revision(git_branch git_commit_hash)

set(COMPILER_FLAGS "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${CMAKE_CXX_FLAGS}")

site_name(BUILD_computerName)

if( BUILD_VF_DOUBLE_ACCURACY )
    set(BUILD_PRECISION "double")
else()
    set(BUILD_PRECISION "float")
endif()

get_target_property(BUILD_COMPILE_OPTIONS project_options INTERFACE_COMPILE_OPTIONS)
get_target_property(BUILD_COMPILE_DEFINITIONS project_options INTERFACE_COMPILE_DEFINITIONS)
get_target_property(BUILD_COMPILE_WARNINGS project_warnings INTERFACE_COMPILE_OPTIONS)

configure_file(${buildInfoInput} ${buildInfoPath}/${buildInfoFile})

set(MY_SRCS ${MY_SRCS} ${buildInfoPath}/${buildInfoFile})
source_group("" FILES  ${buildInfoPath}/${buildInfoFile})
