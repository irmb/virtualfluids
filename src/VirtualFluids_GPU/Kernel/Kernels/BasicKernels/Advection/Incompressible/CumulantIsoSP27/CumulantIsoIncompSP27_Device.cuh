#ifndef LB_KERNEL_CUM_ISO_TEST_INCOMP_SP_27
#define LB_KERNEL_CUM_ISO_TEST_INCOMP_SP_27

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_Cum_IsoTest_Incomp_SP_27(real omega,
	unsigned int* bcMatD,
	unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	real* DDStart,
	real* dxxUx,
	real* dyyUy,
	real* dzzUz,
	int size_Mat,
	bool EvenOrOdd);

#endif