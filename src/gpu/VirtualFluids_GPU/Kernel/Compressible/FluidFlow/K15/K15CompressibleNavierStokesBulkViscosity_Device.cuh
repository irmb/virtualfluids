#ifndef K15CompressibleNavierStokesBulkViscosity_Device_H
#define K15CompressibleNavierStokesBulkViscosity_Device_H

#include <DataTypes.h>
#include <curand.h>


__global__ void K15CompressibleNavierStokesBulkViscosity_Device(real omega,
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