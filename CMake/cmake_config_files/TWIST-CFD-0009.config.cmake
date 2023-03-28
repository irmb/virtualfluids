SET(CMAKE_CUDA_ARCHITECTURES "86")
SET(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)

list(APPEND USER_APPS "apps/gpu/LBM/ActuatorLine")
list(APPEND USER_APPS "apps/gpu/LBM/SphereScaling")