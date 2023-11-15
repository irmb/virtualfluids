#ifndef LB_KERNEL_BGK_INCOMP_SP_27_H
#define LB_KERNEL_BGK_INCOMP_SP_27_H

#include <DataTypes.h>
#include <curand.h>

__global__ void B92IncompressibleNavierStokes_Device(
    real omega,
    unsigned int* bcMatD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* DDStart,
    int size_Mat,
    bool EvenOrOdd);

#endif