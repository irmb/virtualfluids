include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/GMock/Link.cmake)
linkGMock(${targetName})

include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/JsonCpp/Link.cmake)
linkJsonCpp(${targetName})