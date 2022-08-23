#ifndef LB_KERNEL_CUMULANT_K15_BULK_COMP_H
#define LB_KERNEL_CUMULANT_K15_BULK_COMP_H

#include <DataTypes.h>
#include <curand.h>


__global__ void LB_Kernel_CumulantK15BulkComp(real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
	bool EvenOrOdd);

#endif