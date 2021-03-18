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

namespace VF
{
namespace LBM
{

inline __host__ __device__ real getDensity(const real *const &f /*[27]*/)
{
    return ((f[DIR::TNE] + f[DIR::BSW]) + (f[DIR::TSE] + f[DIR::BNW])) + ((f[DIR::BSE] + f[DIR::TNW]) + (f[DIR::TSW] + f[DIR::BNE])) +
           (((f[DIR::NE] + f[DIR::SW]) + (f[DIR::SE] + f[DIR::NW])) + ((f[DIR::TE] + f[DIR::BW]) + (f[DIR::BE] + f[DIR::TW])) +
            ((f[DIR::BN] + f[DIR::TS]) + (f[DIR::TN] + f[DIR::BS]))) +
           ((f[DIR::E] + f[DIR::W]) + (f[DIR::N] + f[DIR::S]) + (f[DIR::T] + f[DIR::B])) + f[DIR::REST];
}


inline __host__ __device__ real getIncompressibleVelocityX1(const real *const &f /*[27]*/)
{
    return ((((f[DIR::TNE] - f[DIR::BSW]) + (f[DIR::TSE] - f[DIR::BNW])) + ((f[DIR::BSE] - f[DIR::TNW]) + (f[DIR::BNE] - f[DIR::TSW]))) +
            (((f[DIR::BE] - f[DIR::TW]) + (f[DIR::TE] - f[DIR::BW])) + ((f[DIR::SE] - f[DIR::NW]) + (f[DIR::NE] - f[DIR::SW]))) + (f[DIR::E] - f[DIR::W]));
}


inline __host__ __device__ real getIncompressibleVelocityX2(const real *const &f /*[27]*/)
{
    return ((((f[DIR::TNE] - f[DIR::BSW]) + (f[DIR::BNW] - f[DIR::TSE])) + ((f[DIR::TNW] - f[DIR::BSE]) + (f[DIR::BNE] - f[DIR::TSW]))) +
            (((f[DIR::BN] - f[DIR::TS]) + (f[DIR::TN] - f[DIR::BS])) + ((f[DIR::NW] - f[DIR::SE]) + (f[DIR::NE] - f[DIR::SW]))) + (f[DIR::N] - f[DIR::S]));
}


inline __host__ __device__ real getIncompressibleVelocityX3(const real *const &f /*[27]*/)
{
    return ((((f[DIR::TNE] - f[DIR::BSW]) + (f[DIR::TSE] - f[DIR::BNW])) + ((f[DIR::TNW] - f[DIR::BSE]) + (f[DIR::TSW] - f[DIR::BNE]))) +
            (((f[DIR::TS] - f[DIR::BN]) + (f[DIR::TN] - f[DIR::BS])) + ((f[DIR::TW] - f[DIR::BE]) + (f[DIR::TE] - f[DIR::BW]))) + (f[DIR::T] - f[DIR::B]));
}



}
}

#endif
