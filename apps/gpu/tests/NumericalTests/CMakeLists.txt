project(NumericalTestsGPU)
set(collectTestFiles ON)
vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluids_GPU GridGenerator basics fftw3 gtest gmock gmock_main cudart_static)

vf_get_library_name(library_name)
groupTarget(${library_name} ${testFolder})

if(BUILD_SHARED_LIBS)
    # add compile option according to
    # https://github.com/google/googletest/blob/master/googletest/README.md#as-a-shared-library-dll
    set_target_properties(${library_name}
            PROPERTIES
            COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")
endif()

linkCUDA()
