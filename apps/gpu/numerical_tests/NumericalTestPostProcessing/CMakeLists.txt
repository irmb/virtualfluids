project(NumericalTestPostProcessing)

vf_add_library(BUILDTYPE binary PRIVATE_LINK gpu_core basics cudart_static)

vf_get_library_name(library_name)
group_target(${library_name} ${testFolder})

target_include_directories(NumericalTestPostProcessing PRIVATE ${CMAKE_CURRENT_LIST_DIR})
