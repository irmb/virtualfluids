#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Cuda/Link.cmake)
#linkCuda(${targetName})
#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/MPI/Link.cmake)
#linkMPI(${targetName})
#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Boost/Link.cmake)
#linkBoost(${targetName} "serialization")


IF(MSVC)
    set(libsToLink ws2_32 GridGenerator basics Traffic) # ws_32 throws an error on Phoenix
ELSE(MSVC)
    #set(libsToLink GridGenerator VirtualFluidsBasics Core Traffic)
    set(libsToLink GridGenerator basics)
ENDIF(MSVC)


vf_add_library(BUILDTYPE shared DEPENDS ${libsToLink})

include (${CMAKE_PATH}/3rd/cuda.cmake)
include (${CMAKE_PATH}/3rd/mpi.cmake)
include (${CMAKE_PATH}/3rd/boost.cmake)
linkBoost ("serialization")

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/gpu/GridGenerator")
target_include_directories(${library_name} PRIVATE "${CMAKE_SOURCE_DIR}/src/basics/Core")



#set(includeDirectories ${CMAKE_SOURCE_DIR}/src/${targetName}
#        ${CMAKE_SOURCE_DIR}/src
#        ${CMAKE_SOURCE_DIR}/src/GridGenerator
#        ${CMAKE_SOURCE_DIR}/src/VirtualFluidsBasics
#        ${CMAKE_SOURCE_DIR}/src/Core
#        ${CMAKE_SOURCE_DIR}/src/Traffic
#        )

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

IF(MSVC)
    set_target_properties(${library_name} PROPERTIES LINK_FLAGS "/ignore:4251")
ENDIF(MSVC)
