
IF(${USE_GCC})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} rt)
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} gomp)
ENDIF()

IF(${USE_INTEL})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} irc)
   SET(CAB_ADDITIONAL_LINK_FLAGS ${CAB_ADDITIONAL_LINK_FLAGS} parallel)
ENDIF()

IF(${USE_ZOLTAN})
   INCLUDE_DIRECTORIES(${ZOLTAN_INCLUDEDIR})
   SET(LINK_LIBRARY optimized ${ZOLTAN_RELEASE_LIBRARY} debug ${ZOLTAN_DEBUG_LIBRARY})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})
ENDIF()

IF(${USE_METIS})
   SET(LINK_LIBRARY optimized ${METIS_RELEASE_LIBRARY} debug ${METIS_DEBUG_LIBRARY})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})
ENDIF()

IF(${USE_VTK})
   SET(LINK_LIBRARY optimized ${VTK_LIBRARIES} debug ${VTK_LIBRARIES})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})
ENDIF()

IF(${USE_CATALYST})
   SET(LINK_LIBRARY optimized vtkPVPythonCatalyst debug vtkPVPythonCatalyst )
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})
   SET(LINK_LIBRARY optimized vtkParallelMPI debug vtkParallelMPI )
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})
ENDIF()


IF(${USE_PYTHON})
   SET(LINK_LIBRARY optimized ${PYTHON_LIBRARY} debug ${PYTHON_LIBRARY})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})
ENDIF()

IF(${USE_FETOL})
   SET(LINK_LIBRARY optimized ${BOND_RELEASE_LIBRARY} debug ${BOND_DEBUG_LIBRARY})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})

   SET(LINK_LIBRARY optimized ${JAVA_JVM_LIBRARY} debug ${JAVA_JVM_LIBRARY})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})

   SET(LINK_LIBRARY optimized ${YAML_RELEASE_LIBRARY} debug ${YAML_DEBUG_LIBRARY})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})

   SET(LINK_LIBRARY optimized ${FETOL_RELEASE_LIBRARY} debug ${FETOL_DEBUG_LIBRARY})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} ${LINK_LIBRARY})
ENDIF()

IF(${USE_DEM_COUPLING})
   INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../DemCoupling/DemCoupling.cmake)
ENDIF()

vf_add_library(BUILDTYPE static DEPENDS basics muparser ${MPI_LIBRARY} ${CAB_ADDITIONAL_LINK_LIBRARIES})

vf_get_library_name(library_name)

target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Connectors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Data)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Interactors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LBM)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Grid)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Visitors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/CoProcessors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Utilities)


target_include_directories(${library_name} PUBLIC ${METIS_INCLUDEDIR})


IF(${USE_BOOST})
   target_include_directories(${library_name} PRIVATE ${Boost_INCLUDE_DIR})
ENDIF()


target_include_directories(${library_name} PRIVATE ${ZOLTAN_INCLUDEDIR})
IF(${USE_VTK})
   target_include_directories(${library_name} PRIVATE ${VTK_INCLUDE_DIRS})
ENDIF()
IF(${USE_FETOL})
   # INCLUDE_DIRECTORIES(${SOURCE_DIR}/VirtualFluidsCore/FETOL)  TODO: Did not exists?
   target_include_directories(${library_name} PRIVATE ${YAML_INCLUDEDIR})
   target_include_directories(${library_name} PRIVATE ${BOND_INCLUDEDIR})
   target_include_directories(${library_name} PRIVATE ${FETOL_INCLUDEDIR})
ENDIF()

