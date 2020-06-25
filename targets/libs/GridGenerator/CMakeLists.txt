setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR})

set(linkDirectories "")
set(libsToLink VirtualFluidsBasics Core)

set(includeDirectories ${CMAKE_SOURCE_DIR}/src/${targetName}
                       ${CMAKE_SOURCE_DIR}/src
                       ${CMAKE_SOURCE_DIR}/src/VirtualFluidsBasics
                       ${CMAKE_SOURCE_DIR}/src/Core )

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${libraryFolder})

set_property(TARGET ${targetName} PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
