#ifndef K17CompressibleNavierStokesWale_Device_H
#define K17CompressibleNavierStokesWale_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K17CompressibleNavierStokesWale_Device(
	real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* neighborWSB,
	real* veloX,
	real* veloY,
	real* veloZ,
	real* DDStart,
	real* turbulentViscosity,
	int size_Mat,
	int level,
	unsigned int timestep, 
	real* forces,
    real* quadricLimiters,
	bool EvenOrOdd);

#endif