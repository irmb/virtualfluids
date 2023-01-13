project(lbm_cuda LANGUAGES CUDA CXX)

vf_add_library(NAME lbm_cuda BUILDTYPE static PUBLIC_LINK basics FOLDER "${CMAKE_SOURCE_DIR}/src/lbm;${CMAKE_SOURCE_DIR}/src/lbm/constants" MODULEFOLDER "lbm")

set_target_properties(lbm_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)

set_source_files_properties(../KernelParameter.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(../CumulantChimera.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(../BGK.cpp PROPERTIES LANGUAGE CUDA)
