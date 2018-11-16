setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR})

set(linkDirectories "")
set(libsToLink GksGpu Core)

set(includeDirectories ${CMAKE_SOURCE_DIR}/src/${targetName}
                       ${CMAKE_SOURCE_DIR}/src
                       ${CMAKE_SOURCE_DIR}/src/GksGpu
                       ${CMAKE_SOURCE_DIR}/src/Core )

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${gksLibraryFolder})

#Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
