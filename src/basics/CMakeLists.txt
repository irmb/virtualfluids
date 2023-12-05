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
#! \author Soeren Peters
#=======================================================================================

include(buildInfo.cmake)

vf_add_library(PUBLIC_LINK logger EXCLUDE buildInfo.in.cpp PRIVATE_LINK project_options project_warnings)

vf_get_library_name (library_name)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry3d)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry3d/KdTree)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/container)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/memory)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/objects)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/transmitter)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utilities)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/writer)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})



IF(MSVC)
    target_compile_definitions(${library_name} PUBLIC NOMINMAX) # Disable Min/Max-Macros
ENDIF(MSVC)

if(VF_ENABLE_BOOST)
    target_link_libraries(${library_name} PRIVATE Boost::boost)
endif()

vf_add_tests()


## WARNINGS
# TODO: Remove this and fix those warnings: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/163
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    target_compile_options(basics PUBLIC "-Wold-style-cast;-Wno-old-style-cast")
    target_compile_options(basics PUBLIC "-Wsign-conversion;-Wno-sign-conversion")
    target_compile_options(basics PUBLIC "-Wdouble-promotion;-Wno-double-promotion")
    target_compile_options(basics PUBLIC "-Wshadow;-Wno-shadow")
    target_compile_options(basics PRIVATE "-Woverloaded-virtual;-Wno-overloaded-virtual")

    target_compile_options(basics PRIVATE "-Wextra") # Wno-unused-parameter comes with -Wextra
    target_compile_options(basics PRIVATE "-Wall") # Wno-unused-function comes with -Wall
    target_compile_options(basics PRIVATE "-Wunused") # -Wno-unused-function comes with -Wunused
    target_compile_options(basics PUBLIC "-Wno-unused-parameter")
    target_compile_options(basics PUBLIC "-Wno-unused-function")
    target_compile_options(basics PRIVATE "-Wno-reorder")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(basics PUBLIC "-Wuseless-cast;-Wno-useless-cast")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(basics PRIVATE "/wd4244") # conversion
    target_compile_options(basics PRIVATE "/wd4100") # unreferenced formal parameter
    target_compile_options(basics PRIVATE "/wd4701") # potentially uninitialized local variable
    target_compile_options(basics PRIVATE "/wd4703") # potentially uninitialized local pointer
    target_compile_options(basics PRIVATE "/wd4458") # declaration of 'XXX' hides class member
    target_compile_options(basics PRIVATE "/wd4305") #'initializing': truncation from double to float
    target_compile_options(basics PRIVATE "/wd4242") #'conversion from 'int' to 'char', possible loss of data
endif()