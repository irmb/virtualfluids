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

namespace VF
{
namespace LBM
{

__host__ __device__ real getDensity(const real *const &f /*[27]*/);

__host__ __device__ real getIncompressibleVelocityX1(const real *const &f /*[27]*/);

__host__ __device__ real getIncompressibleVelocityX2(const real *const &f /*[27]*/);

__host__ __device__ real getIncompressibleVelocityX3(const real *const &f /*[27]*/);

}
}

#endif
