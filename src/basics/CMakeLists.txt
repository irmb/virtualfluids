
include(Core/buildInfo.cmake)

vf_add_library(PUBLIC_LINK logger MPI::MPI_CXX EXCLUDE buildInfo.in.cpp)

vf_get_library_name (library_name)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Core)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry3d)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry3d/KdTree)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/container)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/memory)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/objects)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/transmitter)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utilities)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/writer)



IF(MSVC)
    target_compile_definitions(${library_name} PUBLIC NOMINMAX) # Disable Min/Max-Macros
ENDIF(MSVC)

if(BUILD_USE_BOOST)
    target_link_libraries(${library_name} PRIVATE Boost::boost)
endif()

vf_add_tests()
