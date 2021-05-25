project(GridGenerator LANGUAGES CUDA CXX)


vf_add_library(PRIVATE_LINK basics OpenMP::OpenMP_CXX)

vf_get_library_name(library_name)
set_target_properties(${library_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# we want to suppress all cuda warnings so far for this library.
target_compile_options(${library_name} PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe "-w" >)
