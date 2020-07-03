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

IF(${USE_BOOST})
   INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${METIS_INCLUDEDIR})
INCLUDE_DIRECTORIES(${ZOLTAN_INCLUDEDIR})
IF(${USE_VTK})
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
ENDIF()
IF(${USE_FETOL})
    INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/FETOL)
    INCLUDE_DIRECTORIES(${YAML_INCLUDEDIR})
    INCLUDE_DIRECTORIES(${BOND_INCLUDEDIR})
    INCLUDE_DIRECTORIES(${FETOL_INCLUDEDIR})
ENDIF()
