setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR}) 

set(linkDirectories "")
set(libsToLink ws2_32)
set(includeDirectories ${CMAKE_SOURCE_DIR}/src/${targetName} ${CMAKE_SOURCE_DIR}/src)

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${libraryFolder})
set_target_properties(${targetName} PROPERTIES LINK_FLAGS "/ignore:4251")

#Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
