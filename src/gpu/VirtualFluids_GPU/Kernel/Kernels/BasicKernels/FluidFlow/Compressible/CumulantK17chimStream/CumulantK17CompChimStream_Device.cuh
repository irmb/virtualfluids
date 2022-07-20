#ifndef LB_Kernel_CUMULANT_K17_COMP_CHIM_SPARSE_H
#define LB_Kernel_CUMULANT_K17_COMP_CHIM_SPARSE_H

#include <DataTypes.h>
#include <curand.h>

__global__ void LB_Kernel_CumulantK17CompChimStream(
	real omega,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	unsigned long size_Mat,
	int level,
	real* forces,
	real* quadricLimiters,
	bool isEvenTimestep,
	const uint* fluidNodeIndices,
	uint numberOfFluidNodes);
#endif
