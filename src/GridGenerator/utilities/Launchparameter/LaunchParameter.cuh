#ifndef kernelHelper_CUH
#define kernelHelper_CUH

#include "GridGenerator/global.h"

#include <GridGenerator_EXPORT.h>

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

class LaunchParameter
{
public:
	HOST GridGenerator_EXPORT LaunchParameter();

	HOST GridGenerator_EXPORT static LaunchParameter make_2D1D_launchParameter(int size, int threadDim);
	HOST GridGenerator_EXPORT static LaunchParameter make_1D1D_launchParameter(int size, int threadDim);

	DEVICE static int getGlobalIdx_2D_1D();
	DEVICE static int getGlobalIdx_1D_1D();

	HOST void print() const;

	dim3 threads;
	dim3 blocks;
};


#endif
