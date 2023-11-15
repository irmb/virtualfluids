
vf_add_library(PUBLIC_LINK cpu_core basics muparser lbm)

set_target_properties(simulationconfig PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(simulationconfig PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})