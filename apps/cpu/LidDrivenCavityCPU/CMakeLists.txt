

PROJECT(LidDrivenCavity)

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluidsCore basics FILES LidDrivenCavity.cpp)

vf_get_library_name (library_name)
target_include_directories(${library_name} PRIVATE ${APPS_ROOT_CPU})