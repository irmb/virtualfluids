project(GksGpu LANGUAGES CUDA CXX)

vf_add_library(PRIVATE_LINK basics lbm_cuda GksMeshAdapter OpenMP::OpenMP_CXX MPI::MPI_CXX)

target_include_directories(GksGpu PRIVATE "${VF_THIRD_DIR}/cuda_samples/")

if (NOT MSVC)
    target_compile_options(GksGpu PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-fPIC>")
endif()

