include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Cuda/Link.cmake)
linkCuda(${targetName})
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/MPI/Link.cmake)
linkMPI(${targetName})
#include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Boost/Link.cmake)
#linkBoost(${targetName} "serialization")