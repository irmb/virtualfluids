#ifndef M02IncompressibleNavierStokes_Device_H
#define M02IncompressibleNavierStokes_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void M02IncompressibleNavierStokes_Device(
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	bool EvenOrOdd);

#endif