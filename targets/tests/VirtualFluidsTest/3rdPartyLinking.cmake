include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/GMock/Link.cmake)
linkGMock(${targetName})
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Boost/Link.cmake)
linkBoost(${targetName} "signals;thread;serialization;filesystem")