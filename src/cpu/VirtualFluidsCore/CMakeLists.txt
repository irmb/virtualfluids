
set(VF_LIBRARIES)

IF(${VFCPU_USE_METIS} AND METIS_RELEASE_LIBRARY AND METIS_DEBUG_LIBRARY)
   list(APPEND VF_LIBRARIES optimized ${METIS_RELEASE_LIBRARY} debug ${METIS_DEBUG_LIBRARY})
ELSE()
   list(APPEND VF_LIBRARIES metis)
ENDIF()

IF(${VFCPU_USE_VTK})
   list(APPEND VF_LIBRARIES optimized ${VTK_LIBRARIES} debug ${VTK_LIBRARIES})
ENDIF()

IF(${VFCPU_USE_CATALYST})
   list(APPEND VF_LIBRARIES optimized vtkPVPythonCatalyst debug vtkPVPythonCatalyst)
   list(APPEND VF_LIBRARIES optimized vtkParallelMPI debug vtkParallelMPI )
ENDIF()

IF(${VFCPU_USE_DEM_COUPLING})
   INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../DemCoupling/DemCoupling.cmake)
ENDIF()

if(BUILD_USE_OPENMP)
    list(APPEND VF_LIBRARIES OpenMP::OpenMP_CXX)
endif()

IF(${VFCPU_USE_LIGGGHTS})
   list(APPEND VF_LIBRARIES optimized ${LIGGGHTS_RELEASE_LIBRARY} debug ${LIGGGHTS_DEBUG_LIBRARY})
ENDIF()

vf_add_library(BUILDTYPE static PUBLIC_LINK basics muparser ${VF_LIBRARIES} PRIVATE_LINK lbm mpi logger)

vf_add_tests()

vf_get_library_name(library_name)

target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Connectors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Data)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Interactors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LBM)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Simulation)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Visitors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/SimulationObservers)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Utilities)


IF(${VFCPU_USE_METIS} AND METIS_INCLUDEDIR)
   target_include_directories(${library_name} PUBLIC ${METIS_INCLUDEDIR})
ENDIF()

target_include_directories(${library_name} PRIVATE ${ZOLTAN_INCLUDEDIR})

IF(${VFCPU_USE_VTK})
   target_include_directories(${library_name} PRIVATE ${VTK_INCLUDE_DIRS})
ENDIF()

IF(${VFCPU_ENABLE_LiggghtsCoupling})
   target_include_directories(${library_name} PUBLIC ${LIGGGHTS_SOURCE_DIR})
ENDIF()
