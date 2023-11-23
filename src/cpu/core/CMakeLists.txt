
set(VF_LIBRARIES)

IF(METIS_RELEASE_LIBRARY AND METIS_DEBUG_LIBRARY)
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

if(BUILD_USE_OPENMP)
    list(APPEND VF_LIBRARIES OpenMP::OpenMP_CXX)
endif()

set(library_name cpu_core)

vf_add_library(NAME ${library_name} BUILDTYPE static PUBLIC_LINK basics muparser ${VF_LIBRARIES} parallel PRIVATE_LINK lbm logger)

vf_add_tests(NAME ${library_name})

target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Connectors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Data)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Interactors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LBM)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LBM/Interpolation)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Simulation)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Visitors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/SimulationObservers)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Utilities)


IF(${VFCPU_USE_VTK})
   target_include_directories(${library_name} PRIVATE ${VTK_INCLUDE_DIRS})
ENDIF()

# TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/139
# if(BUILD_USE_MPI)
    target_link_libraries(${library_name} PRIVATE MPI::MPI_CXX)
# endif()
