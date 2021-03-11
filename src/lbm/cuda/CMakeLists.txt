project(lbmCuda LANGUAGES CUDA CXX)


vf_add_library(NAME lbmCuda BUILDTYPE static PUBLIC_LINK basics FOLDER ../../lbm)

#set_target_properties(lbm PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
set_target_properties(lbmCuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(lbm PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_source_files_properties(../CalcMac.cpp PROPERTIES LANGUAGE CUDA)
