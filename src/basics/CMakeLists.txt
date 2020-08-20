
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/geometry3d)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/geometry3d/Creator)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/geometry3d/KdTree)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/container)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/memory)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/objects)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/transmiitter)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/utilities)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/writer)


include (${CMAKE_SOURCE_DIR}/CMake/VirtualFluidsMacros.cmake)
#SET(CAB_ADDITIONAL_LINK_LIBRARIES ${MPI_CXX_LIBRARIES})

include(Core/buildInfo.cmake)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/VirtualFluidsCore)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

IF(MSVC)
    add_compile_definitions(NOMINMAX)                # Disable Min/Max-Macros
ENDIF(MSVC)

vf_add_library(BUILDTYPE static
               FOLDER
        ${CMAKE_CURRENT_LIST_DIR}/geometry3d
        ${CMAKE_CURRENT_LIST_DIR}/geometry3d/creator
        ${CMAKE_CURRENT_LIST_DIR}/geometry3d/KdTree
        ${CMAKE_CURRENT_LIST_DIR}/basics/container
        ${CMAKE_CURRENT_LIST_DIR}/basics/memory
        ${CMAKE_CURRENT_LIST_DIR}/basics/objects
        ${CMAKE_CURRENT_LIST_DIR}/basics/transmiitter
        ${CMAKE_CURRENT_LIST_DIR}/basics/utilities
        ${CMAKE_CURRENT_LIST_DIR}/basics/writer
        ${CMAKE_CURRENT_LIST_DIR}/Core
        ${CMAKE_CURRENT_LIST_DIR}/Core/Logger
        ${CMAKE_CURRENT_LIST_DIR}/Core/Logger/implementations
        ${CMAKE_CURRENT_LIST_DIR}/Core/Input
        ${CMAKE_CURRENT_LIST_DIR}/Core/Input/ConfigData
        ${CMAKE_CURRENT_LIST_DIR}/Core/Input/ConfigFileReader
        ${CMAKE_CURRENT_LIST_DIR}/Core/Input/ConfigInput
        ${CMAKE_CURRENT_LIST_DIR}/Core/Input/JsonInput
        ${CMAKE_CURRENT_LIST_DIR}/Core/StringUtilities
        ${CMAKE_CURRENT_LIST_DIR}/Core/Timer

        EXCLUDE buildInfo.in.cpp
        )

vf_get_library_name (library_name)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/Core)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

include (${CMAKE_PATH}/3rd/mpi.cmake)

vf_add_tests()