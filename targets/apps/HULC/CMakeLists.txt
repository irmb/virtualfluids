setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR}) 

set(linkDirectories "")
set(libsToLink Logger Input GridGenerator VirtualFluids_GPU utilities)
set(includeDirectories "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/src/Logger" "${CMAKE_SOURCE_DIR}/src/Input" "${CMAKE_SOURCE_DIR}/src/GridGenerator" "${CMAKE_SOURCE_DIR}/src/VirtualFluids_GPU" "${CMAKE_SOURCE_DIR}/src/utilities")

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

buildExe(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${appFolder})

# Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)