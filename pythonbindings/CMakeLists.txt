project(VirtualFluidsPython LANGUAGES CUDA CXX)
IF(BUILD_VF_GPU)
    pybind11_add_module(pyfluids src/VirtualFluidsModulesGPU.cpp)
    set_source_files_properties(src/VirtualFluidsModulesGPU.cpp PROPERTIES LANGUAGE CUDA)

    target_link_libraries(pyfluids PRIVATE GridGenerator VirtualFluids_GPU basics lbmCuda logger)
    target_include_directories(pyfluids PRIVATE ${VF_THIRD_DIR}/cuda_samples/)

ENDIF()
IF(BUILD_VF_CPU)
    pybind11_add_module(pyfluids src/VirtualFluidsModulesCPU.cpp)
    pybind11_add_module(pymuparser src/muParser.cpp)

    # TODO: Move this to MuParser CMakeLists.txt
    set_target_properties(muparser PROPERTIES POSITION_INDEPENDENT_CODE ON)

    target_compile_definitions(pyfluids PRIVATE VF_METIS VF_MPI)
    target_compile_definitions(pymuparser PRIVATE VF_METIS VF_MPI)

    target_link_libraries(pyfluids PRIVATE simulationconfig VirtualFluidsCore muparser basics)
    target_link_libraries(pymuparser PRIVATE muparser)
ENDIF()
target_include_directories(pyfluids PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(pyfluids PRIVATE ${CMAKE_BINARY_DIR})