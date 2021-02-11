# Note: CMake support is community-based. The maintainers do not use CMake
# internally.

cmake_minimum_required(VERSION 2.8.12)

if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(googletest-distribution)
set(GOOGLETEST_VERSION 1.10.0)

##########################################################
# Enable MACOSX_RPATH support for gtest and gmock
# to avoid CMP0042 warnings
IF(APPLE)
  SET(CMAKE_MACOSX_RPATH 1)
ENDIF()
##########################################################

if (CMAKE_VERSION VERSION_LESS "3.1")
  add_definitions(-std=c++11)
else()
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  if(NOT CYGWIN)
    set(CMAKE_CXX_EXTENSIONS OFF)
  endif()
endif()

enable_testing()

include(CMakeDependentOption)
include(GNUInstallDirs)

#Note that googlemock target already builds googletest
#option(BUILD_GMOCK "Builds the googlemock subproject" ON)
#option(INSTALL_GTEST "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)" OFF)

#if(BUILD_GMOCK)
  add_subdirectory( googlemock )
#else()
#  add_subdirectory( googletest )
#endif()

groupTarget(gmock ${thirdFolder}/googletest)
groupTarget(gmock_main ${thirdFolder}/googletest)
groupTarget(gtest ${thirdFolder}/googletest)
groupTarget(gtest_main ${thirdFolder}/googletest)