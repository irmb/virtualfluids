project(GridGenerator LANGUAGES CXX)

vf_add_library(PRIVATE_LINK basics OpenMP::OpenMP_CXX)

if(NOT MSVC) 
   target_compile_options(GridGenerator PRIVATE "-Wno-strict-aliasing")
endif()