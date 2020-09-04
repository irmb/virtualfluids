
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/IncludsList.cmake)

IF(${USE_GCC})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} rt)
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} gomp)
ENDIF()

IF(${USE_INTEL})
   SET(CAB_ADDITIONAL_LINK_LIBRARIES ${CAB_ADDITIONAL_LINK_LIBRARIES} irc)
   SET(CAB_ADDITIONAL_LINK_FLAGS ${CAB_ADDITIONAL_LINK_FLAGS} parallel)
ENDIF()


vf_add_library(BUILDTYPE static DEPENDS basics MuParser ${CAB_ADDITIONAL_LINK_LIBRARIES})

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/geometry3d)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/basics/container)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/basics/memory)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/basics/objects)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/basics/transmitter)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/basics/utilities)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/basics/writer)

