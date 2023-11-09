#ifndef LB_KERNEL_AD_INCOMP_7_H
#define LB_KERNEL_AD_INCOMP_7_H

#include <DataTypes.h>
#include <curand.h>

__global__ void LB_Kernel_AD_Incomp_7(real diffusivity,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* DD7,
	unsigned long long numberOfLBnodes,
	bool EvenOrOdd);

#endif