setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR})

set(linkDirectories "")
set(libsToLink "VirtualFluidsBasics" "VirtualFluids")
set(includeDirectories  ${PE_ROOT}/src  ${PE_BINARY_DIR}/src ${CMAKE_SOURCE_DIR}/src)
#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${libraryFolder})
generateExportHeader(${targetName} ${CMAKE_SOURCE_DIR}/src/${targetName})
set_target_properties(${targetName} PROPERTIES LINK_FLAGS "/ignore:4251")

#Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
