#ifndef LB_KERNEL_CUMULANT_K15_SPONGE_COMP_H
#define LB_KERNEL_CUMULANT_K15_SPONGE_COMP_H

#include <DataTypes.h>
#include <curand.h>

__global__ void LB_Kernel_CumulantK15SpongeComp(real omega,
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