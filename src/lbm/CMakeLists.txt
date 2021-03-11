project(lbm LANGUAGES CXX)

if(BUILD_VF_GPU)
    enable_language(CUDA)
endif()

vf_add_library(BUILDTYPE shared PUBLIC_LINK basics)
vf_add_library(NAME lbmStatic BUILDTYPE static PUBLIC_LINK basics)

if(BUILD_VF_GPU)
    #set_target_properties(lbm PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    set_target_properties(lbmStatic PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #set_target_properties(lbm PROPERTIES POSITION_INDEPENDENT_CODE ON)
    set_source_files_properties(CalcMac.cpp PROPERTIES LANGUAGE CUDA)
endif()


vf_add_tests()

