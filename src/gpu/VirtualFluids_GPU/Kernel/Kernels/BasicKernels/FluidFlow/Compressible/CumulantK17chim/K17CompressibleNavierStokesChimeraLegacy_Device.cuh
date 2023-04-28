#ifndef LB_Kernel_CUMULANT_K17_COMP_CHIM_H
#define LB_Kernel_CUMULANT_K17_COMP_CHIM_H

#include <DataTypes.h>
#include <curand.h>

__global__ void K17CompressibleNavierStokesChimeraLegacy_Device(
	real omega,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	int size_Mat,
	int level,
	bool bodyForce,
	real* forces,
	real* bodyForceX,
	real* bodyForceY,
	real* bodyForceZ,
	real* quadricLimiters,
	bool isEvenTimestep);
#endif
