#ifndef C06CompressibleNavierStokes_Device_H
#define C06CompressibleNavierStokes_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void C06CompressibleNavierStokes_Device(
	real s9,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	bool EvenOrOdd);

#endif