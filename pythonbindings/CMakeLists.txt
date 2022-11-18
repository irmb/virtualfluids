project(VirtualFluidsPython LANGUAGES CUDA CXX)
IF(BUILD_VF_GPU)
    pybind11_add_module(python_bindings MODULE src/VirtualFluidsModulesGPU.cpp)
    set_target_properties(python_bindings PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pythonbindings/pyfluids
                        OUTPUT_NAME "bindings")
    set_source_files_properties(src/VirtualFluidsModulesGPU.cpp PROPERTIES LANGUAGE CUDA)

    target_link_libraries(python_bindings PRIVATE GridGenerator VirtualFluids_GPU basics lbmCuda logger)
    target_include_directories(python_bindings PRIVATE ${VF_THIRD_DIR}/cuda_samples/)

ENDIF()
IF(BUILD_VF_CPU)
    pybind11_add_module(python_bindings MODULE src/VirtualFluidsModulesCPU.cpp)
    pybind11_add_module(pymuparser MODULE src/muParser.cpp)

    set_target_properties(python_bindings PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pythonbindings/pyfluids
    OUTPUT_NAME "bindings")

    set_target_properties(pymuparser PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pythonbindings/pyfluids
    OUTPUT_NAME "pymuparser_bindings")

    # TODO: Move this to MuParser CMakeLists.txt
    set_target_properties(muparser PROPERTIES POSITION_INDEPENDENT_CODE ON)

    target_compile_definitions(python_bindings PRIVATE VF_METIS VF_MPI)
    target_compile_definitions(pymuparser PRIVATE VF_METIS VF_MPI)

    target_link_libraries(python_bindings PRIVATE simulationconfig VirtualFluidsCore muparser basics)
    target_link_libraries(pymuparser PRIVATE muparser)
ENDIF()
target_include_directories(python_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(python_bindings PRIVATE ${CMAKE_BINARY_DIR})