#directory pathes for header files

set(VirtualFluidsCore_source_dir ${SOURCE_PATH}/cpu/VirtualFluidsCore)
vf_get_library_name(library_name)

set (SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/cpu/")
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/BoundaryConditions)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/Connectors)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/Data)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/Interactors)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/LBM)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/Parallel)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/Grid)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/Visitors)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/CoProcessors)
INCLUDE_DIRECTORIES(${VirtualFluidsCore_source_dir}/Utilities)

INCLUDE_DIRECTORIES(${THIRD_PATH})


