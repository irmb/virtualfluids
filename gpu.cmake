#############################################################
###                  Core                                 ###
#############################################################

add_subdirectory(src/gpu/GridGenerator)
add_subdirectory(src/gpu/VirtualFluids_GPU)

if(BUILD_VF_ALL_SAMPLES)
    list(APPEND USER_APPS 
    "apps/gpu/LBM/ActuatorLine"
    "apps/gpu/LBM/SphereScaling" 
    "apps/gpu/LBM/TGV_3D")
endif()

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
    if(NOT BUILD_VF_UNIT_TESTS) # in this case googletest is already included.
        add_subdirectory(${VF_THIRD_DIR}/googletest)
    endif()

    add_subdirectory(3rdParty/fftw/fftw-3.3.7)
    add_subdirectory(apps/gpu/tests/NumericalTests)
    add_subdirectory(apps/gpu/tests/NumericalTestPostProcessing)
endif()
