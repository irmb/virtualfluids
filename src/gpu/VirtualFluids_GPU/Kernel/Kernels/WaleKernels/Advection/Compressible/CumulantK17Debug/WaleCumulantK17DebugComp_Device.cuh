#ifndef LB_KERNEL_WALE_CUMULANT_K17_DEBUG_COMP_H
#define LB_KERNEL_WALE_CUMULANT_K17_DEBUG_COMP_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_WaleCumulantK17DebugComp(
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
    real* quadricLimiters,
	bool EvenOrOdd);

#endif