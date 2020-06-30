
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/numerics/geometry3d)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/numerics/geometry3d/Creator)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/numerics/geometry3d/KdTree)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/basics/container)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/basics/memory)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/basics/objects)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/basics/transmiitter)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/basics/utilities)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsBasic/basics/writer)


include (${CMAKE_SOURCE_DIR}/CMake/VirtualFluidsMacros.cmake)
#message (${MPI_LIBRARY})
SET(CAB_ADDITIONAL_LINK_LIBRARIES ${MPI_LIBRARY})

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/VirtualFluidsCore)
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

vf_add_library(BUILDTYPE static
        DEPENDS core
               FOLDER
        ${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d
        ${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d/Creator
        ${CMAKE_CURRENT_LIST_DIR}/numerics/geometry3d/KdTree
        ${CMAKE_CURRENT_LIST_DIR}/basics/container
        ${CMAKE_CURRENT_LIST_DIR}/basics/memory
        ${CMAKE_CURRENT_LIST_DIR}/basics/objects
        ${CMAKE_CURRENT_LIST_DIR}/basics/transmiitter
        ${CMAKE_CURRENT_LIST_DIR}/basics/utilities
        ${CMAKE_CURRENT_LIST_DIR}/basics/writer
        )
vf_get_library_name (library_name)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../core)