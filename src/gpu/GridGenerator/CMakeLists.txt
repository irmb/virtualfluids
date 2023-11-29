project(GridGenerator LANGUAGES CXX)

vf_add_library(PUBLIC_LINK basics logger OpenMP::OpenMP_CXX)
vf_add_tests()

if(NOT MSVC) 
   target_compile_options(GridGenerator PRIVATE "-Wno-strict-aliasing")
endif()

if(BUILD_VF_UNIT_TESTS)
   target_include_directories(GridGeneratorTests PRIVATE "${VF_ROOT_DIR}/src/basics/")
   set_source_files_properties(geometries/BoundingBox/BoundingBoxTest.cpp PROPERTIES LANGUAGE CUDA)
   set_source_files_properties(grid/MultipleGridBuilderFacadeTest.cpp PROPERTIES LANGUAGE CUDA)
endif()