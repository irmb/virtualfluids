project(GksGpu LANGUAGES CUDA CXX)

vf_add_library(PRIVATE_LINK basics lbmCuda GksMeshAdapter OpenMP::OpenMP_CXX MPI::MPI_CXX)

target_include_directories(GksGpu PRIVATE "${VF_THIRD_DIR}/cuda_samples/")
target_compile_options(GksGpu PRIVATE "-fPIC")