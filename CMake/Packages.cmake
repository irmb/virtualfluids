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
include(FetchContent)

# spdlog
set(spdlog_version "v1.9.1")
set(spdlog_url "https://github.com/gabime/spdlog")
message(STATUS "Fetching spdlog: ${spdlog_version}")
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY ${spdlog_url}
    GIT_TAG ${spdlog_version})

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
    find_package(OpenMP REQUIRED)
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
    set(pybind_version "v2.10.4")
    set(pybind_url "https://github.com/pybind/pybind11")
    message(STATUS "Fetching pybind: ${pybind_version}")
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY ${pybind_url}
        GIT_TAG ${pybind_version})

    FetchContent_MakeAvailable(pybind11)
endif()