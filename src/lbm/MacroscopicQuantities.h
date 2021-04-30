#ifndef LBM_CALCMAC_H
#define LBM_CALCMAC_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/Core/DataTypes.h>

#include "D3Q27.h"

namespace vf
{
namespace lbm
{

inline __host__ __device__ real getDensity(const real *const &f /*[27]*/)
{
    return ((f[dir::TNE] + f[dir::BSW]) + (f[dir::TSE] + f[dir::BNW])) + ((f[dir::BSE] + f[dir::TNW]) + (f[dir::TSW] + f[dir::BNE])) +
           (((f[dir::NE] + f[dir::SW]) + (f[dir::SE] + f[dir::NW])) + ((f[dir::TE] + f[dir::BW]) + (f[dir::BE] + f[dir::TW])) +
            ((f[dir::BN] + f[dir::TS]) + (f[dir::TN] + f[dir::BS]))) +
           ((f[dir::E] + f[dir::W]) + (f[dir::N] + f[dir::S]) + (f[dir::T] + f[dir::B])) + f[dir::REST];
}


inline __host__ __device__ real getIncompressibleVelocityX1(const real *const &f /*[27]*/)
{
    return ((((f[dir::TNE] - f[dir::BSW]) + (f[dir::TSE] - f[dir::BNW])) + ((f[dir::BSE] - f[dir::TNW]) + (f[dir::BNE] - f[dir::TSW]))) +
            (((f[dir::BE] - f[dir::TW]) + (f[dir::TE] - f[dir::BW])) + ((f[dir::SE] - f[dir::NW]) + (f[dir::NE] - f[dir::SW]))) + (f[dir::E] - f[dir::W]));
}


inline __host__ __device__ real getIncompressibleVelocityX2(const real *const &f /*[27]*/)
{
    return ((((f[dir::TNE] - f[dir::BSW]) + (f[dir::BNW] - f[dir::TSE])) + ((f[dir::TNW] - f[dir::BSE]) + (f[dir::BNE] - f[dir::TSW]))) +
            (((f[dir::BN] - f[dir::TS]) + (f[dir::TN] - f[dir::BS])) + ((f[dir::NW] - f[dir::SE]) + (f[dir::NE] - f[dir::SW]))) + (f[dir::N] - f[dir::S]));
}


inline __host__ __device__ real getIncompressibleVelocityX3(const real *const &f /*[27]*/)
{
    return ((((f[dir::TNE] - f[dir::BSW]) + (f[dir::TSE] - f[dir::BNW])) + ((f[dir::TNW] - f[dir::BSE]) + (f[dir::TSW] - f[dir::BNE]))) +
            (((f[dir::TS] - f[dir::BN]) + (f[dir::TN] - f[dir::BS])) + ((f[dir::TW] - f[dir::BE]) + (f[dir::TE] - f[dir::BW]))) + (f[dir::T] - f[dir::B]));
}



}
}

#endif
