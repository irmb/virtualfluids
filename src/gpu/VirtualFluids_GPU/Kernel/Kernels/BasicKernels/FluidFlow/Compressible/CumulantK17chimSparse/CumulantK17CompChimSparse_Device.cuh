#ifndef LB_Kernel_CUMULANT_K17_COMP_CHIM_SPARSE_H
#define LB_Kernel_CUMULANT_K17_COMP_CHIM_SPARSE_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_CumulantK17CompChimSparse(
	real omega,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	int size_Mat,
	int level,
	real* forces,
	real* quadricLimiters,
	bool isEvenTimestep,
	const uint* fluidNodeIndices,
	uint numberOfFluidNodes);
#endif
