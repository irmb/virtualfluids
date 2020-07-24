#ifndef LB_KERNEL_WALE_CUMULANT_K15_COMP_H
#define LB_KERNEL_WALE_CUMULANT_K15_COMP_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_WaleCumulantK15Comp(real omega,
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
	bool EvenOrOdd);

#endif