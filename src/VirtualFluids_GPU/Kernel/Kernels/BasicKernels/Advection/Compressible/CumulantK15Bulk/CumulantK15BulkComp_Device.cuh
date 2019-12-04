#ifndef LB_KERNEL_CUM_ONE_COMP_BULK_SP_27_H
#define LB_KERNEL_CUM_ONE_COMP_BULK_SP_27_H

#include <DataTypes.h>
#include <curand.h>


extern "C" __global__ void LB_Kernel_CumulantK15BulkComp(real omega,
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