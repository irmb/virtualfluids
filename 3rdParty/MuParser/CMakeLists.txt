

vf_add_library(BUILDTYPE static)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
