#ifndef K17CompressibleNavierStokesBulkViscosity_Device_H
#define K17CompressibleNavierStokesBulkViscosity_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K17CompressibleNavierStokesBulkViscosity_Device(real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
    real* quadricLimiters,
	bool EvenOrOdd);
#endif