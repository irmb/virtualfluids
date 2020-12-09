cmake_minimum_required(VERSION 3.1)
project(simulationconfig)

set(VFBUILDER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB HEADERS
        ${VFBUILDER_INCLUDE_DIR}/simulationconfig/*.h
        ${VFBUILDER_INCLUDE_DIR}/simulationconfig/boundaryconditions/*.h)

file(GLOB SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_library(simulationconfig STATIC ${SRC} ${HEADERS})

set_target_properties(simulationconfig PROPERTIES CXX_STANDARD 11 POSITION_INDEPENDENT_CODE ON)

target_include_directories(simulationconfig PUBLIC ${VFBUILDER_INCLUDE_DIR})
target_link_libraries(simulationconfig PRIVATE VirtualFluidsCore basics muparser)
target_compile_definitions(simulationconfig PRIVATE VF_METIS VF_MPI)
target_include_directories(simulationconfig PUBLIC ${CMAKE_BINARY_DIR})

