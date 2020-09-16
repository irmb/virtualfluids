# Old VF basics defines. TODO: might be deleted in future.
#LIST(APPEND VF_COMPILER_DEFINITION CAB_BOOST)
#LIST(APPEND VF_COMPILER_DEFINITION BOOST_SIGNALS_NO_DEPRECATION_WARNING)
#LIST(APPEND VF_COMPILER_DEFINITION CAB_RUBY)
#LIST(APPEND VF_COMPILER_DEFINITION mpe=mpilog)
#LIST(APPEND VF_COMPILER_DEFINITION noshlib)
#LIST(APPEND VF_COMPILER_DEFINITION SINGLEPRECISION)


include(Core/buildInfo.cmake)

vf_add_library(BUILDTYPE static EXCLUDE buildInfo.in.cpp VirtualFluidsDefinitions.in.h)

vf_get_library_name (library_name)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Core)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry3d)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry3d/KdTree)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/basics/container)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/basics/memory)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/basics/objects)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/basics/parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/basics/transmitter)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/basics/utilities)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/basics/writer)



IF(MSVC)
    target_compile_definitions(${library_name} PUBLIC NOMINMAX) # Disable Min/Max-Macros
ENDIF(MSVC)


linkMPI()

vf_add_tests()