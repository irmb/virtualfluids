project(VirtualFluidsPython)

pybind11_add_module(pyfluids src/VirtualfluidsModule.cpp)

target_link_libraries(pyfluids PRIVATE simulationconfig VirtualFluidsCore muparser basics)

target_include_directories(pyfluids PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(pyfluids PRIVATE ${CMAKE_BINARY_DIR})