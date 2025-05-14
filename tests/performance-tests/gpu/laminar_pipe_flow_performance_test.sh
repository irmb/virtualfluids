#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
source ./tests/performance-tests/__performance_test_executer.sh

# 1. set reference data directory (must match the folder structure in https://git.rz.tu-bs.de/irmb/virtualfluids-reference-data)
REFERENCE_DATA_DIR=performance_tests/gpu

# 2. set cmake flags for the build of VirtualFluids
CMAKE_FLAGS="--preset=make_gpu -DCMAKE_CUDA_ARCHITECTURES=75"

# 3. define the application to be executed
APPLICATION="./build/bin/LaminarPipeFlowGPU ./apps/gpu/LaminarPipeFlowGPU/laminarpipeflow_performance_test.cfg"

# 4. set the path to the produced data
RESULT_DATA_DIR=output/LaminarPipeFlow

# 5. set the name of the metadata file (.yaml)
META_DATA_NAME=LaminarPipeFlow.yaml


run_performance_test "$REFERENCE_DATA_DIR" "$CMAKE_FLAGS" "$APPLICATION" "$RESULT_DATA_DIR" "$META_DATA_NAME"
