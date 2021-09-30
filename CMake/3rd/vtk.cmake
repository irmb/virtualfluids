
#########################################################################
# VTK_DIR needs to bet set to the VTK build directory in the config file.
#########################################################################
find_package(VTK REQUIRED)
vf_get_library_name(library_name)

target_include_directories(${library_name} PRIVATE ${VTK_INCLUDE_DIRS})

target_link_libraries(${library_name} PRIVATE ${VTK_LIBRARIES})
