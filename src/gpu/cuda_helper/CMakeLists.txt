project(cuda LANGUAGES CUDA CXX)

vf_add_library(PUBLIC_LINK logger PRIVATE_LINK project_options project_warnings)

target_include_directories(cuda_helper SYSTEM PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} "${VF_THIRD_DIR}/cuda_samples/")
