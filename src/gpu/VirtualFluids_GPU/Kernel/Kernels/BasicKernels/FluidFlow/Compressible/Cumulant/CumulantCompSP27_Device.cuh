#ifndef LB_KERNEL_CUM_COMP_SP_27_H
#define LB_KERNEL_CUM_COMP_SP_27_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_Cum_Comp_SP_27(real s9,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	bool EvenOrOdd);

#endif