cmake_minimum_required(VERSION 3.0)
project(METIS)

set(GKLIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/GKlib" CACHE PATH "path to GKlib")

# Configure libmetis library.
if(BUILD_SHARED_LIBS)
  set(METIS_LIBRARY_TYPE SHARED)
else()
  set(METIS_LIBRARY_TYPE STATIC)
endif()

include(${GKLIB_PATH}/GKlibSystem.cmake)


add_subdirectory("libmetis")

target_include_directories(metis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/libmetis)
target_include_directories(metis PRIVATE ${GKLIB_PATH})

target_include_directories(metis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


groupTarget(metis ${thirdFolder})