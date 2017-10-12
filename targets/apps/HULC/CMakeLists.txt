setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR}) 

set(linkDirectories "")
set(libsToLink GridGenerator VirtualFluids_GPU VirtualFluidsBasics)
set(includeDirectories "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/src/GridGenerator" "${CMAKE_SOURCE_DIR}/src/VirtualFluids_GPU")

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildExe(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${appFolder})

# Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)