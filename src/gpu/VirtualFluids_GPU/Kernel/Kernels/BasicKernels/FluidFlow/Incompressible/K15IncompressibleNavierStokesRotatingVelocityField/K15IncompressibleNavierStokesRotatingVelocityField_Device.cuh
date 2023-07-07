#ifndef K15IncompressibleNavierStokesRotatingVelocityField_Device_H
#define K15IncompressibleNavierStokesRotatingVelocityField_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K15IncompressibleNavierStokesRotatingVelocityField_Device(
	real omega,
	real deltaPhi,
	real angularVelocity,
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