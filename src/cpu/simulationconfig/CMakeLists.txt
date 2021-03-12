project(simulationconfig)


vf_add_library(NAME simulationconfig PUBLIC_LINK basics muparser MPI::MPI_CXX ${VF_LIBRARIES} PRIVATE_LINK lbm)

set_target_properties(simulationconfig PROPERTIES CXX_STANDARD 14 POSITION_INDEPENDENT_CODE ON)

target_include_directories(simulationconfig PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
