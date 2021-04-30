#################################################################################
#   _    ___      __              __________      _     __
# | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
# | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
# | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
# |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
#
#################################################################################
#  required cmake versions
#  CMAKE 3.13: target_link_options
#  CMAKE 3.15: CMAKE_MSVC_RUNTIME_LIBRARY
#################################################################################
cmake_minimum_required(VERSION 3.15..3.20 FATAL_ERROR)

project(VirtualFluids CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ".cmake")
set(libraryFolder "libs")
set(testFolder    "tests")
set(appFolder     "apps")
set(thirdFolder   "3rd")

set (VF_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set (VF_THIRD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty)
set (VF_SRC_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/src)
set (VF_ROOT_DIR  ${CMAKE_CURRENT_SOURCE_DIR})

#################################################################################
#  OPTIONS
#################################################################################
option(BUILD_VF_CPU "Build VirtualFluids cpu variant" OFF)
option(BUILD_VF_GPU "Build VirtualFluids gpu variant" OFF)

option(BUILD_USE_OPENMP "Build VirtualFluids with openmp" ON)


# vf gpu
option(BUILD_VF_GPU          "Build VirtualFluids GPU"     ON )
option(BUILD_VF_GKS          "Build VirtualFluids GKS"     OFF )
option(BUILD_VF_TRAFFIC      "Build VirtualFluids Traffic" OFF)
option(BUILD_JSONCPP         "Builds json cpp "            OFF)
option(BUILD_NUMERIC_TESTS   "Build numeric tests"         OFF)

option(BUILD_VF_UNIT_TESTS "Build VirtualFluids unit tests" OFF)
option(BUILD_VF_CLANG_TIDY "Add the clang tidy checks to the targets" OFF)
option(BUILD_VF_INCLUDE_WHAT_YOU_USE "Add IWYU to the targets" OFF)
option(BUILD_VF_CPPCHECK "Add cppcheck to the targets" OFF)
option(BUILD_VF_COVERAGE "Add the -coverage compiler flag." OFF)

option(BUILD_SHARED_LIBS "" OFF)
option(BUILD_WARNINGS_AS_ERRORS "" OFF)

# windows: use multi-threaded dynamically-linked runtime library
if(BUILD_SHARED_LIBS)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

option(BUILD_VF_PYTHON_BINDINGS "" OFF)

option(BUILD_VF_DOUBLE_ACCURACY "Use double accuracy" OFF)


#################################################################################
#  MACROS
#################################################################################
include(CMakePrintHelpers)
include(${VF_CMAKE_DIR}/VirtualFluidsMacros.cmake)

IF( BUILD_VF_DOUBLE_ACCURACY )
    list(APPEND VF_COMPILER_DEFINITION VF_DOUBLE_ACCURACY)
ENDIF()

if(BUILD_VF_GPU)
    include(CheckLanguage)
    check_language(CUDA)

    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR "CUDA Compiler was requested but is not found on the system.")
    endif()

    if(MSVC) 
        # With the MSVC compiler we got this warning: nvcc : The -std=c++14 flag is not supported with the configured host compiler. Flag will be ignored.
        # But we build the c++ code with C++14. Until we have not a solution here, we set the standard to 11 here.
        set(CMAKE_CUDA_STANDARD 11)
    else()
        set(CMAKE_CUDA_STANDARD 14)
    endif()

    set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        message(WARNING "CMAKE_CUDA_ARCHITECTURES was not defined and is set to 30 (CUDA support until 10.1 only).")
        set(CMAKE_CUDA_ARCHITECTURES 30)
    endif()

    message("CUDA Architecture: ${CMAKE_CUDA_ARCHITECTURES}")
endif()

#################################################################################
#  COMMON LIBRARIES
#################################################################################
if(BUILD_VF_UNIT_TESTS)
    add_subdirectory(${VF_THIRD_DIR}/googletest)
    include(GoogleTest)
    enable_testing()
endif()

if(BUILD_USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

find_package(MPI REQUIRED)


add_subdirectory(src/basics)


#################################################################################
#  VIRTUAL FLUIDS CPU / GPU
#################################################################################
if (BUILD_VF_CPU)
    include (cpu.cmake)
endif()
if(BUILD_VF_GPU)
    include (gpu.cmake)
endif()