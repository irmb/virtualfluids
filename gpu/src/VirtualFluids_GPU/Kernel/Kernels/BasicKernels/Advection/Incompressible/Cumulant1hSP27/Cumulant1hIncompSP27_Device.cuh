#ifndef LB_KERNEL_CUM_1H_INCOMP_SP_27_H
#define LB_KERNEL_CUM_1H_INCOMP_SP_27_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_Cum_1h_Incomp_SP_27(real omega,
	real deltaPhi,
	real angularVelocity,
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