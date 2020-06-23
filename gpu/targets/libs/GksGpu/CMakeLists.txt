setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR}) 

set(linkDirectories "")
set(libsToLink GksMeshAdapter Core VirtualFluidsBasics)

set(includeDirectories ${CMAKE_SOURCE_DIR}/src/${targetName}
                       ${CMAKE_SOURCE_DIR}/src
                       ${CMAKE_SOURCE_DIR}/src/GksMeshAdapter
                       ${CMAKE_SOURCE_DIR}/src/Core
                       ${CMAKE_SOURCE_DIR}/src/VirtualFluidsBasics )

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${gksLibraryFolder})

#Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
