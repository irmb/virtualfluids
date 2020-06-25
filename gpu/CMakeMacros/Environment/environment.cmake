unset(BUILD_computerName)
site_name(BUILD_computerName)
MESSAGE(STATUS "computer name: " ${BUILD_computerName})

include(${CMAKE_SOURCE_DIR}/MachineFiles/${BUILD_computerName})