include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/MPI/Link.cmake)
linkMPI(${targetName})
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Boost/Link.cmake)
linkBoost(${targetName} "system;date_time;thread;serialization;chrono;regex")
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Metis/Link.cmake)
linkMetis(${targetName})