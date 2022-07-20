#ifndef LB_Kernel_CUMULANT_K17_COMP_CHIM_REDESIGN_H
#define LB_Kernel_CUMULANT_K17_COMP_CHIM_REDESIGN_H

#include <DataTypes.h>
#include <curand.h>

__global__ void LB_Kernel_CumulantK17CompChimRedesigned(
    real omega,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    real* distributions,
    unsigned long numberOfLBnodes,
    int level,
    real* forces,
    real* quadricLimiters,
    real* rho,
    real* veloX,
    real* veloY,
    real* veloZ,
    bool isEvenTimestep,
    const uint* fluidNodeIndices,
    uint numberOfFluidNodes);
#endif
