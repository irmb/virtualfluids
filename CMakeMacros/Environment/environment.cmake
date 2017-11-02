unset(BUILD_computerName)
site_name(BUILD_computerName)
MESSAGE(STATUS "computer name: " ${BUILD_computerName})

include(${CMAKE_CURRENT_LIST_DIR}/MachineFiles/${BUILD_computerName})