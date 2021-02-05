
function(linkCUDA)

    set(CUDA_CUT_INCLUDE_DIR "${VF_THIRD_DIR}/cuda_samples/")

    vf_get_library_name(library_name)
    target_include_directories(${library_name} PRIVATE ${CUDA_CUT_INCLUDE_DIR})
    target_include_directories(${library_name} PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

    # set the following properties only for specific targets
    # set_property(TARGET ${targetName} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    # set_property(TARGET ${targetName} PROPERTY CUDA_64_BIT_DEVICE_CODE ON)
endfunction()
