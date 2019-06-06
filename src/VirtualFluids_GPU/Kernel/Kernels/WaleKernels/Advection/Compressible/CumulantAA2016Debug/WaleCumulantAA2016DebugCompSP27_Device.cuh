#ifndef LB_KERNEL_WALE_CUM_AA2016_DEBUG_COMP_SP_27_H
#define LB_KERNEL_WALE_CUM_AA2016_DEBUG_COMP_SP_27_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_Wale_Cum_AA2016_Debug_Comp_SP_27(
	real omega_in,
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
	real* gSij,
	real* gSDij,
	real* gDxvx,
	real* gDyvx,
	real* gDzvx,
	real* gDxvy,
	real* gDyvy,
	real* gDzvy,
	real* gDxvz,
	real* gDyvz,
	real* gDzvz,
	int size_Mat,
	int level,
	real* forces,
	bool EvenOrOdd);

#endif