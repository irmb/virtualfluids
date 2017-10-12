include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Metis/Link.cmake)
linkMetis(${targetName})
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/MPI/Link.cmake)
linkMPI(${targetName})
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Boost/Link.cmake)
linkBoost(${targetName} "Serialization")
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Cuda/Link.cmake)
linkCuda(${targetName})