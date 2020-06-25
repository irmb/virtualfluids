#ifndef LB_KERNEL_CUMULANT_K15_IMCOMP_H
#define LB_KERNEL_CUMULANT_K15_IMCOMP_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_CumulantK15Incomp(real s9,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	bool EvenOrOdd);

#endif