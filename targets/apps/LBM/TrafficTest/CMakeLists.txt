setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR}) 

set(linkDirectories "")
set(libsToLink Traffic)
set(includeDirectories "${CMAKE_SOURCE_DIR}/src" 
                       "${CMAKE_SOURCE_DIR}/src/Traffic")

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildExe(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${lbmAppFolder})

# Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)