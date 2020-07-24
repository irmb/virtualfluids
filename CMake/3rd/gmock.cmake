
    if(BUILD_SHARED_LIBS)
        add_definitions(-DGTEST_LINKED_AS_SHARED_LIBRARY)
    endif()


    target_include_directories(${targetName} PRIVATE ${GMOCK_ROOT}/googlemock/include)
    target_include_directories(${targetName} PRIVATE ${GMOCK_ROOT}/googletest/include)

    target_link_directories(${targetName} PRIVATE ${GMOCK_ROOT}/build/lib)

    target_link_libraries(${targetName} PRIVATE gtest gmock gmock_main)
