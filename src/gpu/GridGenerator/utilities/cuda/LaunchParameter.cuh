#ifndef kernelHelper_CUH
#define kernelHelper_CUH

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>

#include "global.h"

class LaunchParameter
{
public:
	CUDA_HOST VIRTUALFLUIDS_GPU_EXPORT LaunchParameter();

	CUDA_HOST VIRTUALFLUIDS_GPU_EXPORT static LaunchParameter make_2D1D_launchParameter(int size, int threadDim);
	CUDA_HOST VIRTUALFLUIDS_GPU_EXPORT static LaunchParameter make_1D1D_launchParameter(int size, int threadDim);

	DEVICE static int getGlobalIdx_2D_1D();
	DEVICE static int getGlobalIdx_1D_1D();

	CUDA_HOST void print() const;

	dim3 threads;
	dim3 blocks;
};


#endif
