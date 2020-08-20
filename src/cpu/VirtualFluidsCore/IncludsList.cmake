#directory pathes for header files

set (SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/cpu/")

INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/BoundaryConditions)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/Connectors)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/Data)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/Interactors)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/LBM)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/Parallel)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/Grid)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/Visitors)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/CoProcessors)
INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/Utilities)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/3rdParty)

