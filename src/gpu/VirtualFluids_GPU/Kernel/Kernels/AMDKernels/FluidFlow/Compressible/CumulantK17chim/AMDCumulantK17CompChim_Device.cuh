#ifndef LB_Kernel_AMD_CUMULANT_K17_COMP_CHIM_H
#define LB_Kernel_AMD_CUMULANT_K17_COMP_CHIM_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_AMDCumulantK17CompChim(
	real omega_in,
	real SGSConstant,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint* neighborWSB,
	real* distributions,
	real* veloX,
    real* veloY,
    real* veloZ,
	unsigned long size_Mat,
	int level,
	bool bodyForce,
	real* forces,
	real* bodyForceX,
	real* bodyForceY,
	real* bodyForceZ,
	real* quadricLimiters,
	bool isEvenTimestep);
#endif
