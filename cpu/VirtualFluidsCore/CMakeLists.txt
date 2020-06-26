
#################################################################
###   PACKAGES						                          ###
#################################################################
SET(WITH_SUBFOLDERS_FOR_SG TRUE)

#INCLUDE(${SOURCE_ROOT}/ThirdParty/MuParser/CMakePackage.txt)

#new VirtualFluids
#sINCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/BoundaryConditions/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Connectors/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Data/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Interactors/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/LBM/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Parallel/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Grid/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Visitors/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/CoProcessors/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Utilities/CMakePackage.txt)

#old VirtualFluids
#INCLUDE(${SOURCE_ROOT}/VirtualFluidsBasic/VirtualFluidsBasic.cmake)

#FETOL
IF(${USE_FETOL})
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/FETOL/CMakePackage.txt)
ENDIF()

INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/IncludsList.cmake)

#SET(CAB_ADDITIONAL_LINK_LIBRARIES muParserLib ${Boost_LIBRARIES} ${MPI_LIBRARY})
SET(CAB_ADDITIONAL_LINK_LIBRARIES ${MPI_LIBRARY})

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
   INCLUDE(${SOURCE_ROOT}/DemCoupling/DemCoupling.cmake)
ENDIF()

message("CAB_ADDITIONAL_LINK_LIBRARIES: " ${CAB_ADDITIONAL_LINK_LIBRARIES})
#################################################################
###   CREATE PROJECT                                          ###
#################################################################
#CREATE_CAB_PROJECT(VirtualFluids STATIC)
include (${PROJECT_SOURCE_DIR}/CMake/VirtualFluidsMacros.cmake)

vf_add_library(BUILDTYPE static DEPENDS VirtualFluidsBasic MuParser)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/numerics/geometry3d)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/numerics/geometry3d/Creator)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/numerics/geometry3d/KdTree)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/basics/container)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/basics/memory)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/basics/objects)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/basics/transmiitter)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/basics/utilities)
target_include_directories(${library_name} PRIVATE ${SOURCE_ROOT}/VirtualFluidsBasic/basics/writer)
