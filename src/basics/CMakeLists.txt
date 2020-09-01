
include(Core/buildInfo.cmake)

vf_add_library(BUILDTYPE static EXCLUDE buildInfo.in.cpp)

vf_get_library_name (library_name)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/Core)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/geometry3d)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/geometry3d/KdTree)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/basics/container)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/basics/memory)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/basics/objects)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/basics/transmitter)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/utilities)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/basics/writer)

IF(MSVC)
    target_compile_definitions(${library_name} PUBLIC NOMINMAX) # Disable Min/Max-Macros
ENDIF(MSVC)

include (${CMAKE_PATH}/3rd/mpi.cmake)

vf_add_tests()