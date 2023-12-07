# #######################################################################################
# ____          ____    __    ______     __________   __      __       __        __
# \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
#  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
#   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
#    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
#     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
#      \    \  |    |   ________________________________________________________________
#       \    \ |    |  |  ______________________________________________________________|
#        \    \|    |  |  |         __          __     __     __     ______      _______
#         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
#          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
#           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
#            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
#
#  This file is part of VirtualFluids. VirtualFluids is free software: you can
#  redistribute it and/or modify it under the terms of the GNU General Public
#  License as published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
# #################################################################################
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
