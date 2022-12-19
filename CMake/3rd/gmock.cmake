#################################################################################
#  Links gmock to the current target.
#  Note: gmock has to be build by the project itself (Located in 3rd).
#################################################################################

function (linkGMOCK)
    vf_get_library_test_name(library_name)
    target_link_libraries(${library_name} PRIVATE GTest::gmock_main)

    if(BUILD_SHARED_LIBS)
        # add compile option according to
        # https://github.com/google/googletest/blob/master/googletest/README.md#as-a-shared-library-dll
        set_target_properties(${library_name}
                PROPERTIES
                COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")
    endif()
endfunction()