setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR}) 

set(linkDirectories "")
set(libsToLink "")
set(includeDirectories "${CMAKE_CURRENT_LIST_DIR}/include")
#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${thirdPartyFolder})
generateExportHeader(${targetName} ${CMAKE_SOURCE_DIR}/3rdParty/${targetName})

#Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
