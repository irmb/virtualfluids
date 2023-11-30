
include(buildInfo.cmake)

vf_add_library(PUBLIC_LINK logger EXCLUDE buildInfo.in.cpp PRIVATE_LINK project_options project_warnings)

vf_get_library_name (library_name)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry3d)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry3d/KdTree)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/container)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/memory)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/objects)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/transmitter)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utilities)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/writer)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})



IF(MSVC)
    target_compile_definitions(${library_name} PUBLIC NOMINMAX) # Disable Min/Max-Macros
ENDIF(MSVC)

if(VF_ENABLE_BOOST)
    target_link_libraries(${library_name} PRIVATE Boost::boost)
endif()

vf_add_tests()


## WARNINGS
# TODO: Remove this and fix those warnings
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    target_compile_options(basics PUBLIC "-Wold-style-cast;-Wno-old-style-cast")
    target_compile_options(basics PUBLIC "-Wsign-conversion;-Wno-sign-conversion")
    target_compile_options(basics PUBLIC "-Wdouble-promotion;-Wno-double-promotion")
    target_compile_options(basics PUBLIC "-Wshadow;-Wno-shadow")
    target_compile_options(basics PRIVATE "-Woverloaded-virtual;-Wno-overloaded-virtual")

    target_compile_options(basics PRIVATE "-Wextra") # Wno-unused-parameter comes with -Wextra
    target_compile_options(basics PRIVATE "-Wall") # Wno-unused-function comes with -Wall
    target_compile_options(basics PRIVATE "-Wunused") # -Wno-unused-function comes with -Wunused
    target_compile_options(basics PUBLIC "-Wno-unused-parameter")
    target_compile_options(basics PUBLIC "-Wno-unused-function")
    # target_compile_options(basics PRIVATE "-Wno-implicit-float-conversion")
    # target_compile_options(basics PRIVATE "-Wno-comment")
    # target_compile_options(basics PRIVATE "-Wno-reorder-ctor")
    # target_compile_options(basics PRIVATE "-Wno-uninitialized")
    target_compile_options(basics PRIVATE "-Wno-reorder")
    # target_compile_options(basics PRIVATE "-Wno-strict-aliasing")
    # target_compile_options(basics PRIVATE "-Wno-format-nonliteral")
    # target_compile_options(basics PRIVATE "-Wno-null-dereference")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(basics PUBLIC "-Wuseless-cast;-Wno-useless-cast")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(basics PRIVATE "/wd4244") # conversion
    target_compile_options(basics PRIVATE "/wd4100") # unreferenced formal parameter
    target_compile_options(basics PRIVATE "/wd4701") # potentially uninitialized local variable
    target_compile_options(basics PRIVATE "/wd4703") # potentially uninitialized local pointer
    target_compile_options(basics PRIVATE "/wd4458") # declaration of 'XXX' hides class member
    target_compile_options(basics PRIVATE "/wd4305") #'initializing': truncation from double to float
    target_compile_options(basics PRIVATE "/wd4242") #'conversion from 'int' to 'char', possible loss of data
endif()