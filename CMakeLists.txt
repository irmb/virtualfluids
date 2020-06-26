cmake_minimum_required(VERSION 3.9..3.17 FATAL_ERROR)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(VirtualFluids CXX)


option(BUILD_VF_CPU "Build VirtualFluids cpu variant" ON)
option(BUILD_VF_GPU "Build VirtualFluids gpu variant" OFF)


if (BUILD_VF_CPU)
    add_subdirectory(cpu)
endif()
if(BUILD_VF_GPU)
    add_subdirectory(gpu)
endif()