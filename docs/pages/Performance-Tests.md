
<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->
# Performance Tests

We conduct performance test of simulation based on the LBM in order to ensure that changes in the source code do not reduce time-to-solution.
They can also be used to check, if a change in the source code improves the simulations' speed. In a performance test we re-run a specific simulation and compare the simulation's performance to a reference performance. The reference performance was determined in a previous run of the same simulation. To measure the performance we use the number of lattice node updates per second (NUPS). with some tolerance the test succeeds. Otherwise the test fails. When the simulation is still fast enough (with some tolerance) the test succeeds. Otherwise the test fails.

# Overview
The performance test use a metadata file in the YAML format that is automatically created when running a simulation. This file contains the performance in NUPS that we compare. In order conduct a meaningful comparison, the tested simulation has to be run on the same hardware as the reference simulation. The utilized hardware for a simulation can be checked in the metadata file as well.

The performance tests are specified in the `tests/performance-tests` directory. Each test is specified in a separate bash file. The file name is the name of the test and needs to end with `\_test.sh`.

Each performance test can be executed manually by running the bash script locally. However, all performance tests are also executed automatically during the CI process (Automatically scheduled or triggered manually). The CI pipeline creates a unique CI-Job for each performance test file.

## Adding a new Test

To add a new test, create a new bash file in the `performance-tests` directory. The file name needs to end with `\_test.sh`. The file name is the name of the test. The file is usually arranged as follows:

```bash
#!/bin/bash
source ./tests/performance-tests/__performance_test_executer.sh

# 1. set reference data directory (must match the folder structure in https://git.rz.tu-bs.de/irmb/virtualfluids-reference-data)
REFERENCE_DATA_DIR=performance_tests/gpu

# 2. set cmake flags for the build of VirtualFluids
CMAKE_FLAGS="--preset=make_gpu_release -DCMAKE_CUDA_ARCHITECTURES=75"

# 3. define the application to be executed
APPLICATION="./build/bin/LaminarPipeFlowGPU ./apps/gpu/LaminarPipeFlowGPU/laminarpipeflow_performance_test.cfg"

# 4. set the path to the produced data
RESULT_DATA_DIR=output/LaminarPipeFlow

# 5. set the name of the metadata file (.yaml)
META_DATA_NAME=LaminarPipeFlow.yaml


run_performance_test "$REFERENCE_DATA_DIR" "$CMAKE_FLAGS" "$APPLICATION" "$RESULT_DATA_DIR" "$META_DATA_NAME"

```

The first line sources the performance test executor script. The next lines are the test-specific configuration. The configuration consists of five parts:

1. The reference data directory. This directory contains the reference data for the test. The directory structure needs to match the directory structure in the [reference data](https://git.rz.tu-bs.de/irmb/virtualfluids-reference-data). The corresponding reference data needs to be uploaded first to the reference data repository.
2. The CMake flags for the build of VirtualFluids. The flags are used to build VirtualFluids for the test.
3. The application to be executed. The application is executed by the performance test executor script.
4. The path to the produced data. The produced data is compared to the reference data.
5. The name of the metadata file. The file name should be the same for both freshly run simulation and the reference simulation

To summarize, adding a new test is usually a two-step process:

1. upload some reference data set to [reference data](https://git.rz.tu-bs.de/irmb/virtualfluids-reference-data) repository:
   - clone the reference data repository
   - create a new folder in the `performance_tests` directory and add your data
   - commit and push the changes
2. create a new test file:
   - copy an existing test file (e.g. laminar_pipe_flow_performance_test.sh) and change the name of the file (still needs to end with `\_test.sh`)
   - adjust points 1-5 in the test file to match your test

Done!
