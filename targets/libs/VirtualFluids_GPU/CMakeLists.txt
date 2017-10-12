setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR}) 

set(linkDirectories "")
set(libsToLink GridGenerator utilities Logger VirtualFluidsBasics ws2_32)
set(includeDirectories ${CMAKE_SOURCE_DIR}/src/${targetName} ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/GridGenerator ${CMAKE_SOURCE_DIR}/src/VirtualFluidsBasics)

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
generateExportHeader(${targetName} ${CMAKE_SOURCE_DIR}/src/${targetName})
groupTarget(${targetName} ${libraryFolder})

#Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
