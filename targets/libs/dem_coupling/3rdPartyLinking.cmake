include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/MPI/Link.cmake)
linkMPI(${targetName})
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Boost/Link.cmake)
linkBoost(${targetName} "signals;thread;serialization;filesystem")
include (${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/Metis/Link.cmake)
linkMetis(${targetName})

target_link_libraries(${targetName}  "pe" "core" "blockforest" "domain_decomposition")
