#!/bin/bash

source ./regression-tests/__regression_test_executer.sh


# 1. set reference data directory (must match the folder structure in https://github.com/irmb/test_data)
REFERENCE_DATA_DIR=regression_tests/cpu/FlowAroundCylinder

# 2. set cmake flags for the build of VirtualFluids
CMAKE_FLAGS="--preset=make_cpu -DCMAKE_BUILD_TYPE=Release"

# 3. define the application to be executed
APPLICATION="./build/bin/FlowAroundCylinder ./apps/cpu/FlowAroundCylinder/FlowAroundCylinder.cfg"

# 4. set the path to the produced data
RESULT_DATA_DIR=output/FlowAroundCylinder


run_regression_test "$REFERENCE_DATA_DIR" "$CMAKE_FLAGS" "$APPLICATION" "$RESULT_DATA_DIR"


# compare forces file manually
file1="$RESULT_DATA_DIR/results/forces.txt"
file2="reference_data/$REFERENCE_DATA_DIR/results/forces.txt"

if cmp -s "$file1" "$file2"; then
    printf 'The file "%s" is the same as "%s"\n' "$file1" "$file2"
else
    printf 'The file "%s" is different from "%s"\n' "$file1" "$file2"
    printf $file1
    cat $file1
    printf $file2
    cat $file2
    exit 1
fi
