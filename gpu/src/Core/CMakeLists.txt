include(buildInfo.cmake)


vf_add_library(BUILDTYPE shared EXCLUDE buildInfo.in.cpp)


include (${CMAKE_PATH}/3rd/mpi.cmake)
include (${CMAKE_PATH}/3rd/boost.cmake)
linkBoost("thread;serialization")