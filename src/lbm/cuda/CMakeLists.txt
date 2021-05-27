project(lbmCuda LANGUAGES CUDA CXX)


vf_add_library(NAME lbmCuda BUILDTYPE static PUBLIC_LINK basics FOLDER ../../lbm)


set_target_properties(lbmCuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

set_source_files_properties(../Distribution27.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(../CumulantChimera.cpp PROPERTIES LANGUAGE CUDA)