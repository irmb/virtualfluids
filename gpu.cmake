#############################################################
###                  Options                              ###
#############################################################

option(VF_GPU_ENABLE_NUMERIC_TESTS "Build numeric tests" OFF)

#############################################################
###                  Libraries                            ###
#############################################################

add_subdirectory(src/gpu/cuda_helper)
add_subdirectory(src/gpu/GridGenerator)
add_subdirectory(src/gpu/core)

#############################################################
###                      Apps                             ###
#############################################################

if(VF_ENABLE_ALL_APPS)
    list(APPEND USER_APPS
    "apps/gpu/DrivenCavityMultiGPU"
    "apps/gpu/AtmosphericBoundaryLayer"
    "apps/gpu/ActuatorLine"
    "apps/gpu/SphereMultiGPU" 
    "apps/gpu/TGV_3D"
    )
endif()

add_subdirectory(apps/gpu/DrivenCavity)
add_subdirectory(apps/gpu/SphereInChannel)

#############################################################
###                   Numeric Tests                       ###
#############################################################

if(VF_GPU_ENABLE_NUMERIC_TESTS)
    if(NOT VF_ENABLE_UNIT_TESTS) # in this case googletest is already included.
        add_subdirectory(${VF_THIRD_DIR}/googletest)
    endif()

    add_subdirectory(3rdParty/fftw/fftw-3.3.7)
    add_subdirectory(apps/gpu/tests/NumericalTests)
    add_subdirectory(apps/gpu/tests/NumericalTestPostProcessing)
endif()
