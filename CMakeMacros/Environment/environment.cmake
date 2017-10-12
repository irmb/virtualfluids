unset(computerName)
site_name(computerName)
MESSAGE(STATUS "computer name: " ${computerName})

include(${CMAKE_CURRENT_LIST_DIR}/MachineFiles/${computerName})