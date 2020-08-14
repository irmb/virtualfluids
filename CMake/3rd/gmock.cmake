
    if(BUILD_SHARED_LIBS)
        add_definitions(-DGTEST_LINKED_AS_SHARED_LIBRARY)
    endif()

    vf_get_library_test_name(library_name)
    target_include_directories(${library_name} PRIVATE ${GMOCK_ROOT}/googlemock/include)
    target_include_directories(${library_name} PRIVATE ${GMOCK_ROOT}/googletest/include)

    target_link_directories(${library_name} PRIVATE ${GMOCK_ROOT}/build/lib)

    target_link_libraries(${library_name} PRIVATE gtest gmock gmock_main)

    #add_compile_options (-std=c++11)
    add_definitions("-std=c++11") #TODO: Really necessary?
