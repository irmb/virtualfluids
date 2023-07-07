#ifndef K15IncompressibleNavierStokes_Device_H
#define K15IncompressibleNavierStokes_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K15IncompressibleNavierStokes_Device(
	real s9,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	bool EvenOrOdd);

#endif