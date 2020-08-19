
set(includeDirectories "${CMAKE_SOURCE_DIR}/src/gpu"
                       "${CMAKE_SOURCE_DIR}/src/basics/Core"
                       "${CMAKE_SOURCE_DIR}/src/gpu/VirtualFluids_GPU"
                       "${CMAKE_SOURCE_DIR}/src/gpu/GridGenerator"
                       "${CMAKE_SOURCE_DIR}/src/basics"
                       "${CMAKE_SOURCE_DIR}/src/gpu/GksMeshAdapter"
                       "${CMAKE_SOURCE_DIR}/src/gpu/GksGpu")



PROJECT(LidDrivenCavity)


vf_add_library(BUILDTYPE binary DEPENDS VirtualFluids_GPU basics GridGenerator GksMeshAdapter GksGpu FILES LidDrivenCavity.cpp)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${includeDirectories} )