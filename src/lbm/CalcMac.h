#ifndef LBM_CALCMAC_H
#define LBM_CALCMAC_H

#include "Core/DataTypes.h"

#ifdef __CUDACC__
#pragma message ( "C Preprocessor got here!" )
#include <cuda_runtime.h>
#else
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#endif 

class LBM
{
public:

__host__ __device__ static real getDensity(const real *const &f /*[27]*/);

__host__ __device__ static real getIncompVelocityX1(const real *const &f /*[27]*/);

__host__ __device__ static real getIncompVelocityX2(const real *const &f /*[27]*/);

__host__ __device__ static real getIncompVelocityX3(const real *const &f /*[27]*/);

};

#endif
