#ifndef LBM_CALCMAC_H
#define LBM_CALCMAC_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#endif

#include "Core/DataTypes.h"

class LBM
{
public:

__host__ __device__ static real getDensity(const real *const &f /*[27]*/);

__host__ __device__ static real getIncompVelocityX1(const real *const &f /*[27]*/);

__host__ __device__ static real getIncompVelocityX2(const real *const &f /*[27]*/);

__host__ __device__ static real getIncompVelocityX3(const real *const &f /*[27]*/);

};

#endif
