#!/bin/bash

#################################
# Driven Cavity Regression Test
#################################

# build VirtualFluids accordingly to our specific test scenario.
# in this case adding -DUSER_APPS="apps/gpu/LBM/DrivenCavity to the cmake command is not necessary, because the DrivenCavity is added to VirtualFluids by default.
mkdir -p build
cmake -B build --preset=make_cpu -DCMAKE_BUILD_TYPE=Release #-DUSER_APPS="apps/cpu/FlowAroundCylinder"
cmake --build build --parallel 8

# execute VirtualFluids
./build/bin/FlowAroundCylinder ./apps/cpu/FlowAroundCylinder/cylinder.cfg 

# set the path to the produced data
PATH_TO_DIR=output/FlowAroundCylinder

# set the path to the reference data.
# `regression-tests/reference_data` is fix `regression_tests/gpu/DrivenCavity_2Levels` must match the structure in https://github.com/irmb/test_data:
PATH_TO_REFERENCE_DIR=reference_data/regression_tests/cpu/FlowAroundCylinder_2023_04

# execute fieldcompare (A more comprehensive manual can be found here https://gitlab.com/dglaeser/fieldcompare)
fieldcompare dir $PATH_TO_DIR $PATH_TO_REFERENCE_DIR --include-files "*.vtu"
