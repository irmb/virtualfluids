#ifndef LB_Kernel_CUMULANT_K17_S_H
#define LB_Kernel_CUMULANT_K17_S_H

#include <DataTypes.h>
#include <curand.h>

template< TurbulenceModel turbulenceModel, bool writeMacroscopicVariables, bool applyBodyForce > __global__ void LB_Kernel_CumulantK17Sponge(
	real omega_in,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	real* rho,
	real* vx,
    real* vy,
    real* vz,
	real* turbulentViscosity,
	real SGSconstant,
	unsigned long numberOfLBnodes,
	int level,
	real* forces,
	real* bodyForceX,
	real* bodyForceY,
	real* bodyForceZ,
	real* coordX,
    real* coordY,
    real* coordZ,
	real* quadricLimiters,
	bool isEvenTimestep,
	const uint *fluidNodeIndices,
    uint numberOfFluidNodes);
#endif
