#ifndef F16CompressibleAdvectionDiffusion_Device_H
#define F16CompressibleAdvectionDiffusion_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void F16CompressibleAdvectionDiffusion_Device(
	real omegaDiffusivity,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	real* distributionsAD,
	unsigned long long numberOfLBnodes,
	real* forces,
	bool isEvenTimestep);

#endif