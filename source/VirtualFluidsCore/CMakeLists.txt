CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(VirtualFluids)
  
#################################################################
###   PACKAGES						                        ###
#################################################################
SET(WITH_SUBFOLDERS_FOR_SG TRUE)
#old VirtualFluids
INCLUDE(${SOURCE_ROOT}/ThirdParty/MuParser/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/ThirdParty/MarchingCubes/CMakePackage.txt)
#INCLUDE(${SOURCE_ROOT}/ThirdParty/Basics/objects/CMakePackage.txt)
# INCLUDE(${SOURCE_ROOT}/ThirdParty/Basics/utilities/CMakePackage.txt)
# #INCLUDE(${SOURCE_ROOT}/ThirdParty/Basics/memory/CMakePackage.txt)
# INCLUDE(${SOURCE_ROOT}/ThirdParty/Basics/container/CMakePackage.txt)
# INCLUDE(${SOURCE_ROOT}/ThirdParty/Basics/writer/CMakePackage.txt)
# #INCLUDE(${SOURCE_ROOT}/ThirdParty/Basics/parallel/CMakePackage.txt)
# INCLUDE(${SOURCE_ROOT}/ThirdParty/Basics/transmitter/CMakePackage.txt)

#new VirtualFluids
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/BoundaryCondition/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Connectors/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Data/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Interactors/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/LBM/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Parallel/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Grid/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Visitors/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/CoProcessors/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Utilities/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/numerics/geometry3d/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/numerics/geometry3d/creator/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/numerics/geometry3d/KdTree/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/basics/objects/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/basics/memory/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/basics/utilities/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/basics/container/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/basics/writer/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/ThirdParty/Library/basics/transmitter/CMakePackage.txt)

#FETOL
IF(${USE_FETOL})
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/FETOL/CMakePackage.txt)
ENDIF()

INCLUDE(${SOURCE_ROOT}/IncludsList.cmake)

SET(CAB_ADDITIONAL_LINK_LIBRARIES muParserLib ${Boost_LIBRARIES} ${MPI_LIBRARY})

IF(${USE_GCC})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} rt)
ENDIF()

IF(${USE_INTEL})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} irc)
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

message("CAB_ADDITIONAL_LINK_LIBRARIES: " ${CAB_ADDITIONAL_LINK_LIBRARIES})
#################################################################
###   CREATE PROJECT                                          ###
#################################################################
CREATE_CAB_PROJECT(VirtualFluids STATIC)
