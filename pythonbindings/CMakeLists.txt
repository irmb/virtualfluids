set(PYFLUIDS_LANGUAGES CXX)

if(BUILD_VF_GPU)
    set(PYFLUIDS_LANGUAGES CUDA CXX)
endif()

project(VirtualFluidsPython LANGUAGES ${PYFLUIDS_LANGUAGES})

pybind11_add_module(python_bindings MODULE src/VirtualFluids.cpp)

set_target_properties(  python_bindings PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pythonbindings/pyfluids
                        OUTPUT_NAME "bindings")

target_link_libraries(python_bindings PRIVATE basics logger mpi)

IF(BUILD_VF_GPU)
    set_source_files_properties(src/VirtualFluids.cpp PROPERTIES LANGUAGE CUDA)

    target_include_directories(python_bindings PRIVATE ${VF_THIRD_DIR}/cuda_samples/)
    target_compile_definitions(python_bindings PRIVATE VF_GPU_PYTHONBINDINGS)

    target_link_libraries(python_bindings PRIVATE GridGenerator VirtualFluids_GPU lbmCuda)
ENDIF()

IF(BUILD_VF_CPU)
    target_compile_definitions(python_bindings PRIVATE VF_METIS VF_MPI VF_CPU_PYTHONBINDINGS)
    target_link_libraries(python_bindings PRIVATE simulationconfig VirtualFluidsCore muparser lbm)

    # include bindings for muparsers
    pybind11_add_module(pymuparser MODULE src/muParser.cpp)

    # TODO: Move this to MuParser CMakeLists.txt
    set_target_properties(muparser PROPERTIES POSITION_INDEPENDENT_CODE ON)

    set_target_properties(  pymuparser PROPERTIES
                            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pythonbindings/pymuparser
                            OUTPUT_NAME "bindings")
    target_compile_definitions(pymuparser PRIVATE VF_METIS VF_MPI)
    target_link_libraries(pymuparser PRIVATE muparser)
ENDIF()


target_include_directories(python_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(python_bindings PRIVATE ${CMAKE_BINARY_DIR})