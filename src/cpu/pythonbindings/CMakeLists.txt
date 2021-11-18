project(VirtualFluidsPython)

pybind11_add_module(pyfluids src/VirtualfluidsModule.cpp)
pybind11_add_module(pymuparser src/muParser.cpp)

# TODO: Move this to MuParser CMakeLists.txt
set_target_properties(muparser PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_compile_definitions(pyfluids PRIVATE VF_METIS VF_MPI)
target_compile_definitions(pymuparser PRIVATE VF_METIS VF_MPI)

target_link_libraries(pyfluids PRIVATE simulationconfig VirtualFluidsCore muparser basics)
target_link_libraries(pymuparser PRIVATE muparser)

target_include_directories(pyfluids PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(pyfluids PRIVATE ${CMAKE_BINARY_DIR})