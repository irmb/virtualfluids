#ifndef LB_Kernel_Cumulant_D3Q27All4_H
#define LB_Kernel_Cumulant_D3Q27All4_H

#include <DataTypes.h>
#include <curand.h>

extern "C" __global__ void LB_Kernel_Cumulant_D3Q27All4(	real omega,
															unsigned int* bcMatD,
															unsigned int* neighborX,
															unsigned int* neighborY,
															unsigned int* neighborZ,
															real* DDStart,
															int size_Mat,
															int level,
															real* forces,
                                                            real* quadricLimiters,
															bool EvenOrOdd);
#endif 