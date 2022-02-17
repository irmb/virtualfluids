#ifndef LB_Kernel_TURBULENT_VISCOSITY_CUMULANT_K17_COMP_CHIM_H
#define LB_Kernel_TURBULENT_VISCOSITY_CUMULANT_K17_COMP_CHIM_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_TurbulentViscosityCumulantK17CompChim(
	real omega_in,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	real* rho,
	real* vx,
    real* vy,
    real* vz,
	real* turbulentViscosity,
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
