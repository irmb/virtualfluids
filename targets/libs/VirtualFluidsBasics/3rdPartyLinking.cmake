include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Boost/Link.cmake)
linkBoost(${targetName} "serialization")
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/MPI/Link.cmake)
linkMPI(${targetName})