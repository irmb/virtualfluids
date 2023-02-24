#!/bin/bash

#################################
# Driven Cavity Regression Test
#################################

# build VirtualFluids accordingly to our specific test scenario.
# in this case adding -DUSER_APPS="apps/gpu/LBM/DrivenCavity to the cmake command is not necessary, because the DrivenCavity is added to VirtualFluids by default.
mkdir -p build
cmake -B build --preset=make_gpu -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=75 #-DUSER_APPS="apps/gpu/LBM/DrivenCavity"
cmake --build build --parallel 8

# execute VirtualFluids
./build/bin/DrivenCavity

# set the path to the produced data
PATH_TO_DIR=output/DrivenCavity

# set the path to the reference data.
# `regression-tests/reference_data` is fix `regression_tests/gpu/DrivenCavity_2Levels` must match the structure in https://github.com/irmb/test_data:
PATH_TO_REFERENCE_DIR=regression-tests/reference_data/regression_tests/gpu/DrivenCavity_2Levels

# execute fieldcompare (A more comprehensive manual can be found here https://gitlab.com/dglaeser/fieldcompare)
fieldcompare dir $PATH_TO_DIR $PATH_TO_REFERENCE_DIR --include-files "*.vtu"
