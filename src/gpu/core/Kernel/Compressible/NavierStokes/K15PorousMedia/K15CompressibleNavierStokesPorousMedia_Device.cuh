#ifndef K15CompressibleNavierStokesPorousMedia_Device_H
#define K15CompressibleNavierStokesPorousMedia_Device_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K15CompressibleNavierStokesPorousMedia_Device(real omega,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	unsigned long long numberOfLBnodes,
	int level,
	real* forces,
	real porosity,
	real darcy,
	real forchheimer,
	unsigned int sizeOfPorousMedia,
	unsigned int* nodeIdsPorousMedia,
	bool EvenOrOdd);

#endif