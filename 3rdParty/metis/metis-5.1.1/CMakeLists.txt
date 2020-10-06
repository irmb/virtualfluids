cmake_minimum_required(VERSION 3.0)
project(METIS C)

set(GKLIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/GKlib")

# Configure libmetis library.
if(BUILD_SHARED_LIBS)
  set(METIS_LIBRARY_TYPE SHARED)
else()
  set(METIS_LIBRARY_TYPE STATIC)
endif()

include(${GKLIB_PATH}/GKlibSystem.cmake)

# METIS' custom options
option(METIS_IDX64 "enable 64 bit ints" OFF)
option(METIS_REAL64 "enable 64 bit floats (i.e., double)" OFF)

if(METIS_IDX64)
  set(METIS_COPTIONS_IDX "-DIDXTYPEWIDTH=64")
else()
  set(METIS_COPTIONS_IDX "-DIDXTYPEWIDTH=32")
endif()

if(METIS_REAL64)
  set(METIS_COPTIONS_REAL "-DREALTYPEWIDTH=64")
else()
  set(METIS_COPTIONS_REAL "-DREALTYPEWIDTH=32")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${METIS_COPTIONS_IDX} ${METIS_COPTIONS_REAL}")


add_subdirectory("libmetis")

target_include_directories(metis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/libmetis)
target_include_directories(metis PRIVATE ${GKLIB_PATH})

target_include_directories(metis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(metis PUBLIC ${METIS_COPTIONS_IDX})
target_compile_definitions(metis PUBLIC ${METIS_COPTIONS_REAL})

groupTarget(metis ${thirdFolder})