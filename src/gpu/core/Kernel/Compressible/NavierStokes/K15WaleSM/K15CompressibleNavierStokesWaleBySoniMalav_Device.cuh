#ifndef K15CompressibleNavierStokesWaleBySoniMalav_Device_H
#define K15CompressibleNavierStokesWaleBySoniMalav_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K15CompressibleNavierStokesWaleBySoniMalav_Device(real omega_in,
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
	real* forces,
	bool EvenOrOdd);


#endif