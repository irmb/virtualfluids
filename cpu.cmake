

add_subdirectory(${VF_THIRD_DIR}/MuParser)

add_subdirectory(src/cpu/VirtualFluidsCore)

set (APPS_ROOT_CPU "${VF_ROOT_DIR}/apps/cpu/")
include(${APPS_ROOT_CPU}/Applications.cmake)