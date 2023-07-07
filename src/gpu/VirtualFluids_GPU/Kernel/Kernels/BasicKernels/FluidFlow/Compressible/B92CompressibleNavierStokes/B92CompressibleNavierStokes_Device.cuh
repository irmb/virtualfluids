#ifndef B92CompressibleNavierStokes_Device_H
#define B92CompressibleNavierStokes_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void B92CompressibleNavierStokes_Device(
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	bool EvenOrOdd);

#endif