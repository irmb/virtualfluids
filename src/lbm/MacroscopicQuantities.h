#ifndef LBM_CALCMAC_H
#define LBM_CALCMAC_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/Core/DataTypes.h>

#include "constants/NumericConstants.h"
#include "constants/D3Q27.h"

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

/*
* Incompressible Macroscopic Quantities
*/
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



/*
* Compressible Macroscopic Quantities
*/
inline __host__ __device__ real getCompressibleVelocityX1(const real *const &f27, const real& rho)
{
    return getIncompressibleVelocityX1(f27) / (rho + constant::c1o1);
}


inline __host__ __device__ real getCompressibleVelocityX2(const real *const &f27, const real& rho)
{
    return getIncompressibleVelocityX2(f27) / (rho + constant::c1o1);
}


inline __host__ __device__ real getCompressibleVelocityX3(const real *const &f27, const real& rho)
{
    return getIncompressibleVelocityX3(f27) / (rho + constant::c1o1);
}

/*
* Pressure
*/
inline __host__ __device__ real getPressure(const real *const &f27, const real& rho, const real& vx, const real& vy, const real& vz)
{
    return (f27[dir::E] + f27[dir::W] + f27[dir::N] + f27[dir::S] + f27[dir::T] + f27[dir::B] + 
    constant::c2o1 * (f27[dir::NE] + f27[dir::SW] + f27[dir::SE] + f27[dir::NW] + f27[dir::TE] + 
                      f27[dir::BW] + f27[dir::BE] + f27[dir::TW] + f27[dir::TN] + f27[dir::BS] + 
                      f27[dir::BN] + f27[dir::TS]) + 
    constant::c3o1 * (f27[dir::TNE] + f27[dir::TSW] + f27[dir::TSE] + f27[dir::TNW] + 
                      f27[dir::BNE] + f27[dir::BSW] + f27[dir::BSE] + f27[dir::BNW]) -
    rho - (vx * vx + vy * vy + vz * vz) * (constant::c1o1 + rho)) * 
    constant::c1o2 + rho; // times zero for incompressible case                 
                          // Attention: op defined directly to op = 1 ; ^^^^(1.0/op-0.5)=0.5
}

}
}

#endif
