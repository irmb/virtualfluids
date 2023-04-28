#ifndef LB_Kernel_CUMULANT_K17_H
#define LB_Kernel_CUMULANT_K17_H

#include <DataTypes.h>
#include <curand.h>

template< TurbulenceModel turbulenceModel, bool writeMacroscopicVariables, bool applyBodyForce > __global__ void K17CompressibleNavierStokes_Device(
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
    unsigned long long numberOfLBnodes,
    int level,
    real* forces,
    real* bodyForceX,
    real* bodyForceY,
    real* bodyForceZ,
    real* quadricLimiters,
    bool isEvenTimestep,
    const uint *fluidNodeIndices,
    uint numberOfFluidNodes);
#endif
