set(PYFLUIDS_LANGUAGES CXX)

if(BUILD_VF_GPU)
    set(PYFLUIDS_LANGUAGES CUDA CXX)
endif()

project(VirtualFluidsPython LANGUAGES ${PYFLUIDS_LANGUAGES})

add_custom_target(python_bindings)

set(PYFLUIDS_DIR ${SKBUILD_PLATLIB_DIR}/pyfluids)
set(PYMUPRASER_DIR ${SKBUILD_PLATLIB_DIR}/pymuparser)


pybind11_add_module(basics_bindings MODULE src/basics/basics.cpp)
set_target_properties(  basics_bindings PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                        OUTPUT_NAME "basics")
target_link_libraries(basics_bindings PRIVATE basics)
target_include_directories(basics_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(basics_bindings PRIVATE ${CMAKE_BINARY_DIR})
add_dependencies(python_bindings basics_bindings)


pybind11_add_module(logger_bindings MODULE src/logger.cpp)
set_target_properties(  logger_bindings PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                        OUTPUT_NAME "logger")
target_link_libraries(logger_bindings PRIVATE logger)
target_include_directories(logger_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(logger_bindings PRIVATE ${CMAKE_BINARY_DIR})
add_dependencies(python_bindings logger_bindings)


pybind11_add_module(lbm_bindings MODULE src/lbm.cpp)
set_target_properties(  lbm_bindings PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                        OUTPUT_NAME "lbm")
target_link_libraries(lbm_bindings PRIVATE lbm)
target_include_directories(lbm_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(lbm_bindings PRIVATE ${CMAKE_BINARY_DIR})
add_dependencies(python_bindings lbm_bindings)

pybind11_add_module(communicator_bindings MODULE src/communicator.cpp)
set_target_properties(  communicator_bindings PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                        OUTPUT_NAME "communicator")
target_link_libraries(communicator_bindings PRIVATE parallel)
target_include_directories(communicator_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(communicator_bindings PRIVATE ${CMAKE_BINARY_DIR})
target_compile_definitions(communicator_bindings PRIVATE VF_MPI)
add_dependencies(python_bindings communicator_bindings)


IF(BUILD_VF_GPU)
    pybind11_add_module(gpu_bindings MODULE src/gpu/gpu.cpp)
    set_target_properties(  gpu_bindings PROPERTIES
                            LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                            OUTPUT_NAME "gpu")
    set_source_files_properties(src/gpu/gpu.cpp PROPERTIES LANGUAGE CUDA)

    target_include_directories(gpu_bindings PRIVATE ${VF_THIRD_DIR}/cuda_samples/)

    target_link_libraries(gpu_bindings PRIVATE basics GridGenerator VirtualFluids_GPU parallel)

    target_include_directories(gpu_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
    target_include_directories(gpu_bindings PRIVATE ${CMAKE_BINARY_DIR})
    add_dependencies(python_bindings gpu_bindings)

ENDIF()

IF(BUILD_VF_CPU)
    pybind11_add_module(cpu_bindings MODULE src/cpu/cpu.cpp)
    set_target_properties(  cpu_bindings PROPERTIES
                            LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                            OUTPUT_NAME "cpu")
    target_link_libraries(cpu_bindings PRIVATE simulationconfig VirtualFluidsCore muparser)

    target_include_directories(cpu_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
    target_include_directories(cpu_bindings PRIVATE ${CMAKE_BINARY_DIR})

    target_compile_definitions(cpu_bindings PRIVATE VF_METIS VF_MPI)
    add_dependencies(python_bindings cpu_bindings)

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

if(BUILD_VF_DOUBLE_ACCURACY)
IF(BUILD_VF_CPU)
    target_compile_definitions(cpu_bindings PRIVATE VF_DOUBLE_ACCURACY)
endif()
    target_compile_definitions(basics_bindings PRIVATE VF_DOUBLE_ACCURACY)
    target_compile_definitions(logger_bindings PRIVATE VF_DOUBLE_ACCURACY)
    target_compile_definitions(lbm_bindings PRIVATE VF_DOUBLE_ACCURACY)
    IF(BUILD_VF_GPU)
    target_compile_definitions(gpu_bindings PRIVATE VF_DOUBLE_ACCURACY)
    endif()
    target_compile_definitions(communicator_bindings PRIVATE VF_DOUBLE_ACCURACY)
endif()