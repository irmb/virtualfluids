
set(includeDirectories "${CMAKE_SOURCE_DIR}/src" 
                       "${CMAKE_SOURCE_DIR}/src/Traffic"
                       "${CMAKE_SOURCE_DIR}/src/GridGenerator"
                       "${CMAKE_SOURCE_DIR}/src/Core")



vf_add_library(BUILDTYPE binary DEPENDS Traffic)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu")