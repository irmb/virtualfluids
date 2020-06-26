#directory pathes for header files

INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/BoundaryConditions)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/Connectors)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/Data)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/Interactors)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/LBM)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/Parallel)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/Grid)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/Visitors)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/CoProcessors)
INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/Utilities)

INCLUDE_DIRECTORIES(${SOURCE_ROOT}/ThirdParty)

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
    INCLUDE_DIRECTORIES(${SOURCE_ROOT}/VirtualFluidsCore/FETOL)
    INCLUDE_DIRECTORIES(${YAML_INCLUDEDIR})
    INCLUDE_DIRECTORIES(${BOND_INCLUDEDIR})
    INCLUDE_DIRECTORIES(${FETOL_INCLUDEDIR})
ENDIF()
