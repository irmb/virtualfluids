PROJECT(GridGeneratorTest)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator FILES gridGenerator.cpp)

set_source_files_properties(gridGenerator.cpp PROPERTIES LANGUAGE CUDA)

linkCUDA()
