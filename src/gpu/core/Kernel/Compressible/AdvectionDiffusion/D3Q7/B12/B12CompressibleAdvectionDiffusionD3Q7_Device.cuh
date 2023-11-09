#ifndef B12CompressibleAdvectionDiffusionD3Q7_Device_H
#define B12CompressibleAdvectionDiffusionD3Q7_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void B12CompressibleAdvectionDiffusionD3Q7_Device(real diffusivity,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* DD7,
	int size_Mat,
	bool EvenOrOdd);

#endif