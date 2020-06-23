CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

########################################################
## C++ PROJECT                                       ###
########################################################
PROJECT(VirtualFluids)
  
#################################################################
###   PACKAGES						                          ###
#################################################################
SET(WITH_SUBFOLDERS_FOR_SG TRUE)

INCLUDE(${SOURCE_ROOT}/ThirdParty/MuParser/CMakePackage.txt)

#new VirtualFluids
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/BoundaryConditions/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Connectors/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Data/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Interactors/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/LBM/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Parallel/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Grid/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Visitors/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/CoProcessors/CMakePackage.txt)
INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/Utilities/CMakePackage.txt)

#old VirtualFluids
INCLUDE(${SOURCE_ROOT}/VirtualFluidsBasics/VirtualFluidsBasics.cmake)

INCLUDE(${SOURCE_ROOT}/VirtualFluidsCore/IncludsList.cmake)

IF(${USE_GCC})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} rt)
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} gomp)
ENDIF()

IF(${USE_INTEL})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} irc)
   SET(CAB_ADDITIONAL_LINK_FLAGS ${CAB_ADDITIONAL_LINK_FLAGS} parallel)
ENDIF()

message("CAB_ADDITIONAL_LINK_LIBRARIES: " ${CAB_ADDITIONAL_LINK_LIBRARIES})
#################################################################
###   CREATE PROJECT                                          ###
#################################################################
CREATE_CAB_PROJECT(VirtualFluids STATIC)
