PROJECT(GridGeneratorTest)

vf_add_library(BUILDTYPE binary DEPENDS basics VirtualFluids_GPU GridGenerator FILES gridGenerator.cpp)

linkCUDA()
