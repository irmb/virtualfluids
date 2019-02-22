setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR})

set(linkDirectories "")
set(libsToLink "")

set(includeDirectories ${CMAKE_SOURCE_DIR}/src/${targetName}
                       ${CMAKE_SOURCE_DIR}/src )

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${libraryFolder})

# Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
