PROJECT(GridGeneratorTest)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator FILES gridGenerator.cpp)

linkCUDA()
