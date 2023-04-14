# Regression Tests

## Overview
The regression tests are specified in the `regression-tests` directory. Each test is specified in a separate bash file. The file name is the name of the test and needs to end with `_test.sh`.

Each regression test can be executed manually by executing it. However, all regression tests are also executed automatically during the CI process (Automatically scheduled or triggered manually). The CI pipeline creates a unique CI-Job for each regression test file.

## Adding a new Test
To add a new test, create a new bash file in the `regression-tests` directory. The file name needs to end with `_test.sh`. The file name is the name of the test. The file is usually arranged as follows:

```bash
#!/bin/bash
source ./regression-tests/__regression_test_executer.sh

# 1. set reference data directory (must match the folder structure in https://github.com/irmb/test_data)
REFERENCE_DATA_DIR=regression_tests/gpu/DrivenCavity_2Levels

# 2. set cmake flags for the build of VirtualFluids
CMAKE_FLAGS="--preset=make_gpu -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=75"

# 3. define the application to be executed
APPLICATION=./build/bin/DrivenCavity

# 4. set the path to the produced data
PATH_TO_DIR=output/DrivenCavity


run_regression_test "$REFERENCE_DATA_DIR" "$CMAKE_FLAGS" "$APPLICATION" "$PATH_TO_DIR"
```

The first line is the shebang line and specifies the interpreter to be used. The second line sources the regression test executer script. The next lines are the test specific configuration. The configuration consists of four parts:
1. The reference data directory. This directory contains the reference data for the test. The directory structure needs to match the directory structure in the [test_data](https://github.com/irmb/test_data). The corresponding reference data needs to be uploaded first to the test_data repository.
2. The cmake flags for the build of VirtualFluids. The flags are used to build VirtualFluids for the test.
3. The application to be executed. The application is executed by the regression test executer script.
4. The path to the produced data. The produced data is compared to the reference data.


To summarize, adding a new test is usually a two step process:
1. upload some reference data set to [test_data](https://github.com/irmb/test_data) repository
    - clone the test_data repository
    - create a new folder in the `regression_tests` directory and add your data
    - commit and push the changes
2. create a new test file 
    - copy an existing test file (e.g. driven_cavity_test.sh) and change the name of the file (still needs to end with `_test.sh`)
    - adjust the points 1-4 in the test file to match your test

Done!
