#############################################################
###                  Libraries                            ###
#############################################################

add_subdirectory(src/gpu/cuda_helper)
add_subdirectory(src/gpu/GridGenerator)
add_subdirectory(src/gpu/core)

#############################################################
###                      Apps                             ###
#############################################################

if(BUILD_VF_ALL_SAMPLES)
    list(APPEND USER_APPS
    "apps/gpu/DrivenCavityMultiGPU"
    "apps/gpu/ActuatorLine"
    "apps/gpu/SphereScaling" 
    "apps/gpu/TGV_3D"
    "apps/gpu/gridGeneratorTest")
endif()

add_subdirectory(apps/gpu/DrivenCavity)
add_subdirectory(apps/gpu/SphereInChannel)
add_subdirectory(apps/gpu/BoundaryLayer)

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
