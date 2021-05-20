#ifndef LB_Kernel_CUMULANT_K17_UNIFIED_H
#define LB_Kernel_CUMULANT_K17_UNIFIED_H

#include <DataTypes.h>
#include <cuda_runtime.h>


namespace vf
{
namespace gpu 
{

extern "C" __global__ void LB_Kernel_CumulantK17Unified(real omega, unsigned int *bcMatD, unsigned int *neighborX,
                                                        unsigned int *neighborY, unsigned int *neighborZ, real *DDStart,
                                                        int size_Mat, int level, real *forces, real *quadricLimiters,
                                                        bool EvenOrOdd);


}
}

#endif
