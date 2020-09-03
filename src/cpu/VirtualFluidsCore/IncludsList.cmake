set(VirtualFluidsCore_source_dir ${VF_SRC_DIR}/cpu/VirtualFluidsCore)
vf_get_library_name(library_name)

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

INCLUDE_DIRECTORIES(${VF_THIRD_DIR})

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
    # INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/FETOL)  TODO: Did not exists?
    INCLUDE_DIRECTORIES(${YAML_INCLUDEDIR})
    INCLUDE_DIRECTORIES(${BOND_INCLUDEDIR})
    INCLUDE_DIRECTORIES(${FETOL_INCLUDEDIR})
ENDIF()

