project(GridGenerator LANGUAGES CUDA CXX)


vf_add_library(PRIVATE_LINK basics OpenMP::OpenMP_CXX)

vf_get_library_name(library_name)
set_target_properties(${library_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# according to linker error when building static libraries.
# https://stackoverflow.com/questions/50033435/cmake-cuda-separate-compilation-static-lib-link-error-on-windows-but-not-on-ubun
if (NOT BUILD_SHARED_LIBRARY)
    set_target_properties(${library_name} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

# supress warning 3123, 3126: calling a __host__ function from __host__ __device__ is not allowed
# supress warning 3152: __host__ function redeclared with __host__ __device__, hence treated as a __host__ __device__ function
# supress warning 3125, 3127: calling a constexpr __host__ function from a __host__ __device__ function
# supress warning 3311: unrecognized #pragma in device code
target_compile_options(${library_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe "--diag_suppress=3123 --diag_suppress=3126 --diag_suppress=3152 --diag_suppress=3125 --diag_suppress=3127 --diag_suppress=3311">) 
