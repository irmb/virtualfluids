#!/bin/bash
source ./regression-tests/__regression_test_executer.sh


# 1. set reference data directory (must match the folder structure in https://github.com/irmb/test_data)
REFERENCE_DATA_DIR=regression_tests/cpu/GyroidsRow

# 2. set cmake flags for the build of VirtualFluids
CMAKE_FLAGS="--preset=make_cpu -DBUILD_USE_BOOST=ON -DCMAKE_BUILD_TYPE=Release"

# 3. define the application to be executed
APPLICATION="mpiexec -np 8 --allow-run-as-root ./build/bin/GyroidsRow ./apps/cpu/GyroidsRow/GyroidsRow_Re_425_u_001_N_100_compressible_regression.cfg"

# 4. set the path to the produced data
RESULT_DATA_DIR=output/GyroidsRow

apt install libboost-all-dev -y

run_regression_test "$REFERENCE_DATA_DIR" "$CMAKE_FLAGS" "$APPLICATION" "$RESULT_DATA_DIR"
