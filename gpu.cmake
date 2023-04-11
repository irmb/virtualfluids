#############################################################
###                  Core                                 ###
#############################################################

add_subdirectory(src/gpu/GridGenerator)
add_subdirectory(src/gpu/VirtualFluids_GPU)

#############################################################
###                  Apps                                 ###
#############################################################

add_subdirectory(apps/gpu/LBM/DrivenCavity)
add_subdirectory(apps/gpu/LBM/SphereGPU)
add_subdirectory(apps/gpu/LBM/BoundaryLayer)

#############################################################
###                   Numeric Tests                       ###
#############################################################

if(BUILD_NUMERIC_TESTS)

    # PATH_NUMERICAL_TESTS can be passed to cmake e.g. cmake .. -DPATH_NUMERICAL_TESTS=/data/
    if(PATH_NUMERICAL_TESTS)
        LIST(APPEND VF_COMPILER_DEFINITION "PATH_NUMERICAL_TESTS=${PATH_NUMERICAL_TESTS}")
    endif()

    if(NOT BUILD_VF_UNIT_TESTS) # in this case googletest is already included.
        add_subdirectory(${VF_THIRD_DIR}/googletest)
    endif()

    add_subdirectory(3rdParty/fftw/fftw-3.3.7)
    add_subdirectory(apps/gpu/tests/NumericalTests)
    add_subdirectory(apps/gpu/tests/NumericalTestPostProcessing)
endif()

#############################################################
###                 Annas Traffic Sim                     ###
#############################################################

if(BUILD_VF_TRAFFIC)
    add_subdirectory(src/gpu/Traffic)
    add_subdirectory(apps/gpu/LBM/TrafficTest)
endif()
