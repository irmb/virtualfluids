setTargetNameToFolderName(${CMAKE_CURRENT_LIST_DIR}) 

set(linkDirectories "")

IF(MSVC)
  set(libsToLink ws2_32 GridGenerator VirtualFluidsBasics Core Traffic) # ws_32 throws an error on Phoenix
ELSE(MSVC)
  set(libsToLink GridGenerator VirtualFluidsBasics Core Traffic)
ENDIF(MSVC)

set(includeDirectories ${CMAKE_SOURCE_DIR}/src/${targetName}
                       ${CMAKE_SOURCE_DIR}/src
                       ${CMAKE_SOURCE_DIR}/src/GridGenerator
                       ${CMAKE_SOURCE_DIR}/src/VirtualFluidsBasics
                       ${CMAKE_SOURCE_DIR}/src/Core 
					   ${CMAKE_SOURCE_DIR}/src/Traffic 
					   )

#glob files and save in MY_SRCS
include(CMakePackage.cmake)

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

buildLib(${targetName} "${MY_SRCS}" "${linkDirectories}" "${libsToLink}" "${includeDirectories}")
groupTarget(${targetName} ${libraryFolder})

IF(MSVC)
  set_target_properties(${targetName} PROPERTIES LINK_FLAGS "/ignore:4251")
ENDIF(MSVC)

#Specify the linking to 3rdParty libs
include(3rdPartyLinking.cmake)
