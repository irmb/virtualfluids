#ifndef K18CompressibleNavierStokes_Device_H
#define K18CompressibleNavierStokes_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K18CompressibleNavierStokes_Device(	
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* F3,
	int size_Mat,
	int level,
	real* forces,
    real* quadricLimiters,
	bool EvenOrOdd);

#endif 