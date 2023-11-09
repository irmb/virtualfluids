#ifndef LBM_CALCMAC_H
#define LBM_CALCMAC_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <gpu/core/LBM/LB.h>
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include "constants/D3Q27.h"


namespace vf::lbm
{

////////////////////////////////////////////////////////////////////////////////////
//! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa  2015.05.001 ]</b></a>
//!
inline __host__ __device__ real getDensity(const real *const &f /*[27]*/)
{
    return ((((f[dir::dPPP] + f[dir::dMMM]) + (f[dir::dMPM] + f[dir::dPMP])) +
             ((f[dir::dMPP] + f[dir::dPMM]) + (f[dir::dMMP] + f[dir::dPPM]))) +
            (((f[dir::d0MP] + f[dir::d0PM]) + (f[dir::d0MM] + f[dir::d0PP])) +
             ((f[dir::dM0P] + f[dir::dP0M]) + (f[dir::dM0M] + f[dir::dP0P])) +
             ((f[dir::dMP0] + f[dir::dPM0]) + (f[dir::dMM0] + f[dir::dPP0]))) +
            f[dir::d000]) +
           ((f[dir::dM00] + f[dir::dP00]) + (f[dir::d0M0] + f[dir::d0P0]) + (f[dir::d00M] + f[dir::d00P]));
}

/*
* Incompressible Macroscopic Quantities
*/
inline __host__ __device__ real getIncompressibleVelocityX1(const real *const &f /*[27]*/)
{
    return ((((f[dir::dPPP] - f[dir::dMMM]) + (f[dir::dPMP] - f[dir::dMPM])) + ((f[dir::dPMM] - f[dir::dMPP]) + (f[dir::dPPM] - f[dir::dMMP]))) +
            (((f[dir::dP0M] - f[dir::dM0P]) + (f[dir::dP0P] - f[dir::dM0M])) + ((f[dir::dPM0] - f[dir::dMP0]) + (f[dir::dPP0] - f[dir::dMM0]))) + (f[dir::dP00] - f[dir::dM00]));
}

inline __host__ __device__ real getIncompressibleVelocityX2(const real *const &f /*[27]*/)
{
    return ((((f[dir::dPPP] - f[dir::dMMM]) + (f[dir::dMPM] - f[dir::dPMP])) + ((f[dir::dMPP] - f[dir::dPMM]) + (f[dir::dPPM] - f[dir::dMMP]))) +
            (((f[dir::d0PM] - f[dir::d0MP]) + (f[dir::d0PP] - f[dir::d0MM])) + ((f[dir::dMP0] - f[dir::dPM0]) + (f[dir::dPP0] - f[dir::dMM0]))) + (f[dir::d0P0] - f[dir::d0M0]));
}

inline __host__ __device__ real getIncompressibleVelocityX3(const real *const &f /*[27]*/)
{
    return ((((f[dir::dPPP] - f[dir::dMMM]) + (f[dir::dPMP] - f[dir::dMPM])) + ((f[dir::dMPP] - f[dir::dPMM]) + (f[dir::dMMP] - f[dir::dPPM]))) +
            (((f[dir::d0MP] - f[dir::d0PM]) + (f[dir::d0PP] - f[dir::d0MM])) + ((f[dir::dM0P] - f[dir::dP0M]) + (f[dir::dP0P] - f[dir::dM0M]))) + (f[dir::d00P] - f[dir::d00M]));
}

inline __host__ __device__ void getIncompressibleMacroscopicValues(const real *const &f /*[27]*/, real &rho, real &vx1, real &vx2, real &vx3)
{
    rho = getDensity(f);
    vx1 = getIncompressibleVelocityX1(f);
    vx2 = getIncompressibleVelocityX2(f);
    vx3 = getIncompressibleVelocityX3(f);
}



/*
* Compressible Macroscopic Quantities
*/
inline __host__ __device__ real getCompressibleVelocityX1(const real *const &f27, const real& rho)
{
    return getIncompressibleVelocityX1(f27) / (rho + basics::constant::c1o1);
}


inline __host__ __device__ real getCompressibleVelocityX2(const real *const &f27, const real& rho)
{
    return getIncompressibleVelocityX2(f27) / (rho + basics::constant::c1o1);
}


inline __host__ __device__ real getCompressibleVelocityX3(const real *const &f27, const real& rho)
{
    return getIncompressibleVelocityX3(f27) / (rho + basics::constant::c1o1);
}

inline __host__ __device__ void getCompressibleMacroscopicValues(const real *const &f /*[27]*/, real &drho, real& oneOverRho, real &vx1, real &vx2, real &vx3)
{
    drho = getDensity(f);
    vx1 = getIncompressibleVelocityX1(f);
    vx2 = getIncompressibleVelocityX2(f);
    vx3 = getIncompressibleVelocityX3(f);
    oneOverRho = vf::basics::constant::c1o1 / (drho + vf::basics::constant::c1o1);
    vx1 *= oneOverRho;
    vx2 *= oneOverRho;
    vx3 *= oneOverRho;
}

inline __host__ __device__ void getCompressibleMacroscopicValues(const real *const &f /*[27]*/, real &drho, real &vx1, real &vx2, real &vx3)
{
    real oneOverRho;
    getCompressibleMacroscopicValues(f, drho, oneOverRho, vx1, vx2, vx3);
}

inline __host__ __device__ void getCompressibleMacroscopicValues(const Distributions27 &dist /*[27] for all nodes*/, uint nodeIndex, real &drho, real &vx1, real &vx2, real &vx3)
{
    drho = ((dist.f[dir::dPPP][nodeIndex] + dist.f[dir::dMMM][nodeIndex]) +
            (dist.f[dir::dPMP][nodeIndex] + dist.f[dir::dMPM][nodeIndex])) +
           ((dist.f[dir::dPMM][nodeIndex] + dist.f[dir::dMPP][nodeIndex]) +
            (dist.f[dir::dMMP][nodeIndex] + dist.f[dir::dPPM][nodeIndex])) +
           (((dist.f[dir::dPP0][nodeIndex] + dist.f[dir::dMM0][nodeIndex]) +
             (dist.f[dir::dPM0][nodeIndex] + dist.f[dir::dMP0][nodeIndex])) +
            ((dist.f[dir::dP0P][nodeIndex] + dist.f[dir::dM0M][nodeIndex]) +
             (dist.f[dir::dP0M][nodeIndex] + dist.f[dir::dM0P][nodeIndex])) +
            ((dist.f[dir::d0PM][nodeIndex] + dist.f[dir::d0MP][nodeIndex]) +
             (dist.f[dir::d0PP][nodeIndex] + dist.f[dir::d0MM][nodeIndex]))) +
           ((dist.f[dir::dP00][nodeIndex] + dist.f[dir::dM00][nodeIndex]) +
            (dist.f[dir::d0P0][nodeIndex] + dist.f[dir::d0M0][nodeIndex]) +
            (dist.f[dir::d00P][nodeIndex] + dist.f[dir::d00M][nodeIndex])) +
           dist.f[dir::d000][nodeIndex];
    real oneOverRho = vf::basics::constant::c1o1 / (drho + vf::basics::constant::c1o1);

    vx1 = ((((dist.f[dir::dPPP][nodeIndex] - dist.f[dir::dMMM][nodeIndex]) +
             (dist.f[dir::dPMP][nodeIndex] - dist.f[dir::dMPM][nodeIndex])) +
            ((dist.f[dir::dPMM][nodeIndex] - dist.f[dir::dMPP][nodeIndex]) +
             (dist.f[dir::dPPM][nodeIndex] - dist.f[dir::dMMP][nodeIndex]))) +
           (((dist.f[dir::dP0M][nodeIndex] - dist.f[dir::dM0P][nodeIndex]) +
             (dist.f[dir::dP0P][nodeIndex] - dist.f[dir::dM0M][nodeIndex])) +
            ((dist.f[dir::dPM0][nodeIndex] - dist.f[dir::dMP0][nodeIndex]) +
             (dist.f[dir::dPP0][nodeIndex] - dist.f[dir::dMM0][nodeIndex]))) +
           (dist.f[dir::dP00][nodeIndex] - dist.f[dir::dM00][nodeIndex]));
    vx1 *= oneOverRho;

    vx2 = ((((dist.f[dir::dPPP][nodeIndex] - dist.f[dir::dMMM][nodeIndex]) +
             (dist.f[dir::dMPM][nodeIndex] - dist.f[dir::dPMP][nodeIndex])) +
            ((dist.f[dir::dMPP][nodeIndex] - dist.f[dir::dPMM][nodeIndex]) +
             (dist.f[dir::dPPM][nodeIndex] - dist.f[dir::dMMP][nodeIndex]))) +
           (((dist.f[dir::d0PM][nodeIndex] - dist.f[dir::d0MP][nodeIndex]) +
             (dist.f[dir::d0PP][nodeIndex] - dist.f[dir::d0MM][nodeIndex])) +
            ((dist.f[dir::dMP0][nodeIndex] - dist.f[dir::dPM0][nodeIndex]) +
             (dist.f[dir::dPP0][nodeIndex] - dist.f[dir::dMM0][nodeIndex]))) +
           (dist.f[dir::d0P0][nodeIndex] - dist.f[dir::d0M0][nodeIndex]));
    vx2 *= oneOverRho;

    vx3 = ((((dist.f[dir::dPPP][nodeIndex] - dist.f[dir::dMMM][nodeIndex]) +
             (dist.f[dir::dPMP][nodeIndex] - dist.f[dir::dMPM][nodeIndex])) +
            ((dist.f[dir::dMPP][nodeIndex] - dist.f[dir::dPMM][nodeIndex]) +
             (dist.f[dir::dMMP][nodeIndex] - dist.f[dir::dPPM][nodeIndex]))) +
           (((dist.f[dir::d0MP][nodeIndex] - dist.f[dir::d0PM][nodeIndex]) +
             (dist.f[dir::d0PP][nodeIndex] - dist.f[dir::d0MM][nodeIndex])) +
            ((dist.f[dir::dM0P][nodeIndex] - dist.f[dir::dP0M][nodeIndex]) +
             (dist.f[dir::dP0P][nodeIndex] - dist.f[dir::dM0M][nodeIndex]))) +
           (dist.f[dir::d00P][nodeIndex] - dist.f[dir::d00M][nodeIndex]));
    vx3 *= oneOverRho;
}

/*
* Pressure
*/
inline __host__ __device__ real getPressure(const real *const &f27, const real& rho, const real& vx, const real& vy, const real& vz)
{
    return (f27[dir::dP00] + f27[dir::dM00] + f27[dir::d0P0] + f27[dir::d0M0] + f27[dir::d00P] + f27[dir::d00M] + 
    basics::constant::c2o1 * (f27[dir::dPP0] + f27[dir::dMM0] + f27[dir::dPM0] + f27[dir::dMP0] + f27[dir::dP0P] + 
                      f27[dir::dM0M] + f27[dir::dP0M] + f27[dir::dM0P] + f27[dir::d0PP] + f27[dir::d0MM] + 
                      f27[dir::d0PM] + f27[dir::d0MP]) + 
    basics::constant::c3o1 * (f27[dir::dPPP] + f27[dir::dMMP] + f27[dir::dPMP] + f27[dir::dMPP] + 
                      f27[dir::dPPM] + f27[dir::dMMM] + f27[dir::dPMM] + f27[dir::dMPM]) -
    rho - (vx * vx + vy * vy + vz * vz) * (basics::constant::c1o1 + rho)) * 
    basics::constant::c1o2 + rho; // times zero for incompressible case                 
                          // Attention: op defined directly to op = 1 ; ^^^^(1.0/op-0.5)=0.5
}


}

#endif
