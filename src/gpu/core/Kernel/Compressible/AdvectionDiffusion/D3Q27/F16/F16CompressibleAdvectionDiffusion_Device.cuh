#ifndef F16CompressibleAdvectionDiffusion_Device_H
#define F16CompressibleAdvectionDiffusion_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void F16CompressibleAdvectionDiffusion_Device(real diffusivity,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* DD27,
	int size_Mat,
	bool EvenOrOdd);

#endif