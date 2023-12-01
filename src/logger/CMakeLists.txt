
vf_add_library(NAME logger PUBLIC_LINK spdlog PRIVATE_LINK project_warnings)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(_flag "/wd4459") # declaration of uint hides global declaration
    target_compile_options(logger PUBLIC "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_flag}>")
    target_compile_options(logger PUBLIC "$<$<COMPILE_LANGUAGE:CXX>:${_flag}>")
endif()
