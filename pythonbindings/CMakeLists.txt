set(PYFLUIDS_LANGUAGES CXX)

if(VF_ENABLE_GPU)
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


pybind11_add_module(parallel_bindings MODULE src/parallel.cpp)
set_target_properties(  parallel_bindings PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                        OUTPUT_NAME "parallel")

target_link_libraries(parallel_bindings PRIVATE parallel)
target_include_directories(parallel_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_include_directories(parallel_bindings PRIVATE ${CMAKE_BINARY_DIR})
add_dependencies(python_bindings parallel_bindings)

target_link_libraries(basics_bindings PRIVATE project_options)
target_link_libraries(logger_bindings PRIVATE project_options)
target_link_libraries(lbm_bindings PRIVATE project_options)
target_link_libraries(parallel_bindings PRIVATE project_options)

IF(VF_ENABLE_GPU)
    pybind11_add_module(gpu_bindings MODULE src/gpu/gpu.cpp)
    set_target_properties(  gpu_bindings PROPERTIES
                            LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                            OUTPUT_NAME "gpu")
    set_source_files_properties(src/gpu/gpu.cpp PROPERTIES LANGUAGE CUDA)

    target_include_directories(gpu_bindings PRIVATE ${VF_THIRD_DIR}/cuda_samples/)

    target_link_libraries(gpu_bindings PRIVATE basics GridGenerator gpu_core parallel)

    target_include_directories(gpu_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
    target_include_directories(gpu_bindings PRIVATE ${CMAKE_BINARY_DIR})
    add_dependencies(python_bindings gpu_bindings)
    target_link_libraries(gpu_bindings PRIVATE project_options)
ENDIF()

IF(VF_ENABLE_CPU)
    pybind11_add_module(cpu_bindings MODULE src/cpu/cpu.cpp)
    set_target_properties(  cpu_bindings PROPERTIES
                            LIBRARY_OUTPUT_DIRECTORY ${PYFLUIDS_DIR}
                            OUTPUT_NAME "cpu")
    target_link_libraries(cpu_bindings PRIVATE simulationconfig cpu_core muparser)

    target_include_directories(cpu_bindings PRIVATE ${CMAKE_SOURCE_DIR}/src/)
    target_include_directories(cpu_bindings PRIVATE ${CMAKE_BINARY_DIR})

    add_dependencies(python_bindings cpu_bindings)

    # include bindings for muparsers
    pybind11_add_module(pymuparser MODULE src/muParser.cpp)

    # TODO: Move this to MuParser CMakeLists.txt
    set_target_properties(muparser PROPERTIES POSITION_INDEPENDENT_CODE ON)

    set_target_properties(  pymuparser PROPERTIES
                            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/pythonbindings/pymuparser
                            OUTPUT_NAME "bindings")
    target_link_libraries(pymuparser PRIVATE muparser)

    target_link_libraries(pymuparser PRIVATE project_options)
    target_link_libraries(cpu_bindings PRIVATE project_options)
ENDIF()
