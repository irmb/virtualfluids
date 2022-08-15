#!/bin/bash


mkdir -p build
cmake -B build --preset=gpu_make -DCMAKE_CUDA_ARCHITECTURES=75 #-DUSER_APPS="apps/gpu/LBM/DrivenCavity"
cd build && make -j 8 && cd ..


./build/bin/DrivenCavity


PATH_TO_DIR=output/DrivenCavity
PATH_TO_REFERENCE_DIR=regression-tests/reference_data/regression_tests_gpu/DrivenCavity


fieldcompare dir $PATH_TO_DIR --reference $PATH_TO_REFERENCE_DIR --include-files "*.vtu"