#ifndef K15CompressibleNavierStokes_Device_H
#define K15CompressibleNavierStokes_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K15CompressibleNavierStokes_Device(	
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
	bool EvenOrOdd);
#endif