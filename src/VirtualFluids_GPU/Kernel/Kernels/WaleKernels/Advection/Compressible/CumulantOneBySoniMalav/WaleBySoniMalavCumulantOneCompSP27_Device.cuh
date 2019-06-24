#ifndef LB_KERNEL_WALE_BY_SONI_MALAV_CUM_ONE_COMP_SP_27_H
#define LB_KERNEL_WALE_BY_SONI_MALAV_CUM_ONE_COMP_SP_27_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_WaleBySoniMalav_Cum_One_Comp_SP_27(real omega_in,
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
	int size_Mat,
	int level,
	real* forces,
	bool EvenOrOdd);


#endif