#ifndef LB_KERNEL_AD_INCOMP_27_H
#define LB_KERNEL_AD_INCOMP_27_H

#include <DataTypes.h>
#include <curand.h>

__global__ void LB_Kernel_AD_Incomp_27(real diffusivity,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* DD27,
	int size_Mat,
	bool EvenOrOdd);

#endif