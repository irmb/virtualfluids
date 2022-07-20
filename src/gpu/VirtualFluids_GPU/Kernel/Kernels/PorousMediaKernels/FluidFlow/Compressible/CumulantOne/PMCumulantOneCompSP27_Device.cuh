#ifndef LB_KERNEL_PM_CUM_ONE_COMP_SP_27_H
#define LB_KERNEL_PM_CUM_ONE_COMP_SP_27_H

#include <DataTypes.h>
#include <curand.h>

__global__ void LB_Kernel_PM_Cum_One_Comp_SP_27(real omega,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	int size_Mat,
	int level,
	real* forces,
	real porosity,
	real darcy,
	real forchheimer,
	unsigned int sizeOfPorousMedia,
	unsigned int* nodeIdsPorousMedia,
	bool EvenOrOdd);

#endif