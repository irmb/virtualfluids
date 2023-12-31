#!/bin/bash
source ./regression-tests/__regression_test_executer.sh


# 1. set reference data directory (must match the folder structure in https://github.com/irmb/test_data)
REFERENCE_DATA_DIR=regression_tests/gpu/ActuatorLine

# 2. set cmake flags for the build of VirtualFluids
CMAKE_FLAGS="--preset=make_gpu -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=75 -DUSER_APPS=apps/gpu/ActuatorLine"

# 3. define the application to be executed
APPLICATION="./build/bin/ActuatorLine ./apps/gpu/ActuatorLine/actuatorline_regressiontest.cfg"

# 4. set the path to the produced data
RESULT_DATA_DIR=output/ActuatorLine


run_regression_test "$REFERENCE_DATA_DIR" "$CMAKE_FLAGS" "$APPLICATION" "$RESULT_DATA_DIR" --ignore-missing-reference-files

