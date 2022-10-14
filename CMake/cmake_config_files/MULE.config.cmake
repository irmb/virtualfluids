SET(CMAKE_CUDA_ARCHITECTURES "75")

set(GPU_APP "apps/gpu/LBM/")
list(APPEND USER_APPS 
    "${GPU_APP}ActuatorLine"
    "${GPU_APP}BoundaryLayer"
)