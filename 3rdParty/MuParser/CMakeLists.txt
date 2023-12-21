# SPDX-License-Identifier: BSD-2-Clause
# SPDX-FileCopyrightText: Copyright (C) 2004 - 2020 Ingo Berg

# By default, build in Release mode. Must appear before project() command
#if (NOT DEFINED CMAKE_BUILD_TYPE)
#    set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
#endif ()

project(muParserProject)

# Bump versions on release
set(MUPARSER_VERSION_MAJOR 2)
set(MUPARSER_VERSION_MINOR 3)
set(MUPARSER_VERSION_PATCH 2)
set(MUPARSER_VERSION ${MUPARSER_VERSION_MAJOR}.${MUPARSER_VERSION_MINOR}.${MUPARSER_VERSION_PATCH})

# Build options
option(ENABLE_SAMPLES "Build the samples" OFF)
#option(ENABLE_OPENMP "Enable OpenMP for multithreading" ON)
#option(BUILD_SHARED_LIBS "Build shared/static libs" ON)


# Credit: https://stackoverflow.com/questions/2368811/how-to-set-warning-level-in-cmake/3818084
if(MSVC)
    # Force to always compile with W4
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4251 /wd4310 /wd4267") # disable all muparser warnings
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    # Update if necessary
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Wno-unknown-pragmas")
endif()

add_library(muparser
        src/muParserBase.cpp
        src/muParserBytecode.cpp
        src/muParserCallback.cpp
        src/muParser.cpp
        src/muParserDLL.cpp
        src/muParserError.cpp
        src/muParserInt.cpp
        src/muParserTest.cpp
        src/muParserTokenReader.cpp
        )

if(VF_ENABLE_OPENMP)
   target_link_libraries(muparser PUBLIC OpenMP::OpenMP_CXX)
endif()
target_include_directories(muparser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# this compiles the "DLL" interface (C API)
target_compile_definitions(muparser PRIVATE MUPARSER_DLL)

if (BUILD_SHARED_LIBS)
    target_compile_definitions(muparser PRIVATE MUPARSERLIB_EXPORTS)
else ()
    target_compile_definitions(muparser PUBLIC MUPARSER_STATIC)
endif()

if (CMAKE_BUILD_TYPE STREQUAL Debug)
    target_compile_definitions(muparser PRIVATE _DEBUG)
endif ()

if(VF_ENABLE_OPENMP)
    target_compile_definitions(muparser PRIVATE MUP_USE_OPENMP)
endif()
set_target_properties(muparser PROPERTIES
        VERSION ${MUPARSER_VERSION}
        SOVERSION ${MUPARSER_VERSION_MAJOR}
        )

if(ENABLE_SAMPLES)
    add_executable(example1 samples/example1/example1.cpp)
    target_link_libraries(example1 muparser)

    add_executable(example2 samples/example2/example2.c)
    target_link_libraries(example2 muparser)
endif()

#include(CTest)
#enable_testing()

#add_executable (t_ParserTest test/t_ParserTest.cpp)
#target_link_libraries(t_ParserTest muparser)
#add_test (NAME ParserTest COMMAND t_ParserTest)

group_target(muparser ${thirdFolder})