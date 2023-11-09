#ifndef K15CompressibleNavierStokesSponge_Device_H
#define K15CompressibleNavierStokesSponge_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K15CompressibleNavierStokesSponge_Device(real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* coordX,
	real* coordY,
	real* coordZ,
	real* DDStart,
	int size_Mat,
	bool EvenOrOdd);

#endif