
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d/Creator)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d/KdTree)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/container)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/memory)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/objects)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/transmiitter)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/utilities)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/basics/writer)


include (${CMAKE_SOURCE_DIR}/CMake/VirtualFluidsMacros.cmake)
SET(CAB_ADDITIONAL_LINK_LIBRARIES ${MPI_CXX_LIBRARIES})

include(Core/buildInfo.cmake)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/VirtualFluidsCore)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

vf_add_library(BUILDTYPE static
               FOLDER
        ${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d
        ${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d/creator
        ${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d/KdTree
        ${CMAKE_CURRENT_LIST_DIR}/basics/container
        ${CMAKE_CURRENT_LIST_DIR}/basics/memory
        ${CMAKE_CURRENT_LIST_DIR}/basics/objects
        ${CMAKE_CURRENT_LIST_DIR}/basics/transmiitter
        ${CMAKE_CURRENT_LIST_DIR}/basics/utilities
        ${CMAKE_CURRENT_LIST_DIR}/basics/writer
        ${CMAKE_CURRENT_LIST_DIR}/Core

        EXCLUDE buildInfo.in.cpp
        )
vf_get_library_name (library_name)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/Core)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
