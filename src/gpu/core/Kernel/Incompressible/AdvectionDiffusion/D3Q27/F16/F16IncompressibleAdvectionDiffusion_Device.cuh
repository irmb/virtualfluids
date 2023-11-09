#ifndef F16IncompressibleAdvectionDiffusion_Device_H
#define F16IncompressibleAdvectionDiffusion_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void F16IncompressibleAdvectionDiffusion_Device(real diffusivity,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* DD27,
	unsigned long long numberOfLBnodes,
	bool EvenOrOdd);

#endif