#ifndef LBM_CALCMAC_H
#define LBM_CALCMAC_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/DataTypes.h>

#include "constants/NumericConstants.h"
#include "constants/D3Q27.h"


namespace vf::lbm
{

////////////////////////////////////////////////////////////////////////////////////
//! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa  2015.05.001 ]</b></a>
//!

inline __host__ __device__ real getDensity(const real *const &f /*[27]*/)
{
    return ((f[dir::DIR_PPP] + f[dir::DIR_MMM]) + (f[dir::DIR_PMP] + f[dir::DIR_MPM])) + ((f[dir::DIR_PMM] + f[dir::DIR_MPP]) + (f[dir::DIR_MMP] + f[dir::DIR_PPM])) +
           (((f[dir::DIR_PP0] + f[dir::DIR_MM0]) + (f[dir::DIR_PM0] + f[dir::DIR_MP0])) + ((f[dir::DIR_P0P] + f[dir::DIR_M0M]) + (f[dir::DIR_P0M] + f[dir::DIR_M0P])) +
            ((f[dir::DIR_0PM] + f[dir::DIR_0MP]) + (f[dir::DIR_0PP] + f[dir::DIR_0MM]))) +
           ((f[dir::DIR_P00] + f[dir::DIR_M00]) + (f[dir::DIR_0P0] + f[dir::DIR_0M0]) + (f[dir::DIR_00P] + f[dir::DIR_00M])) + f[dir::DIR_000];
}

/*
* Incompressible Macroscopic Quantities
*/
inline __host__ __device__ real getIncompressibleVelocityX1(const real *const &f /*[27]*/)
{
    return ((((f[dir::DIR_PPP] - f[dir::DIR_MMM]) + (f[dir::DIR_PMP] - f[dir::DIR_MPM])) + ((f[dir::DIR_PMM] - f[dir::DIR_MPP]) + (f[dir::DIR_PPM] - f[dir::DIR_MMP]))) +
            (((f[dir::DIR_P0M] - f[dir::DIR_M0P]) + (f[dir::DIR_P0P] - f[dir::DIR_M0M])) + ((f[dir::DIR_PM0] - f[dir::DIR_MP0]) + (f[dir::DIR_PP0] - f[dir::DIR_MM0]))) + (f[dir::DIR_P00] - f[dir::DIR_M00]));
}


inline __host__ __device__ real getIncompressibleVelocityX2(const real *const &f /*[27]*/)
{
    return ((((f[dir::DIR_PPP] - f[dir::DIR_MMM]) + (f[dir::DIR_MPM] - f[dir::DIR_PMP])) + ((f[dir::DIR_MPP] - f[dir::DIR_PMM]) + (f[dir::DIR_PPM] - f[dir::DIR_MMP]))) +
            (((f[dir::DIR_0PM] - f[dir::DIR_0MP]) + (f[dir::DIR_0PP] - f[dir::DIR_0MM])) + ((f[dir::DIR_MP0] - f[dir::DIR_PM0]) + (f[dir::DIR_PP0] - f[dir::DIR_MM0]))) + (f[dir::DIR_0P0] - f[dir::DIR_0M0]));
}


inline __host__ __device__ real getIncompressibleVelocityX3(const real *const &f /*[27]*/)
{
    return ((((f[dir::DIR_PPP] - f[dir::DIR_MMM]) + (f[dir::DIR_PMP] - f[dir::DIR_MPM])) + ((f[dir::DIR_MPP] - f[dir::DIR_PMM]) + (f[dir::DIR_MMP] - f[dir::DIR_PPM]))) +
            (((f[dir::DIR_0MP] - f[dir::DIR_0PM]) + (f[dir::DIR_0PP] - f[dir::DIR_0MM])) + ((f[dir::DIR_M0P] - f[dir::DIR_P0M]) + (f[dir::DIR_P0P] - f[dir::DIR_M0M]))) + (f[dir::DIR_00P] - f[dir::DIR_00M]));
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

/*
* Pressure
*/
inline __host__ __device__ real getPressure(const real *const &f27, const real& rho, const real& vx, const real& vy, const real& vz)
{
    return (f27[dir::DIR_P00] + f27[dir::DIR_M00] + f27[dir::DIR_0P0] + f27[dir::DIR_0M0] + f27[dir::DIR_00P] + f27[dir::DIR_00M] + 
    basics::constant::c2o1 * (f27[dir::DIR_PP0] + f27[dir::DIR_MM0] + f27[dir::DIR_PM0] + f27[dir::DIR_MP0] + f27[dir::DIR_P0P] + 
                      f27[dir::DIR_M0M] + f27[dir::DIR_P0M] + f27[dir::DIR_M0P] + f27[dir::DIR_0PP] + f27[dir::DIR_0MM] + 
                      f27[dir::DIR_0PM] + f27[dir::DIR_0MP]) + 
    basics::constant::c3o1 * (f27[dir::DIR_PPP] + f27[dir::DIR_MMP] + f27[dir::DIR_PMP] + f27[dir::DIR_MPP] + 
                      f27[dir::DIR_PPM] + f27[dir::DIR_MMM] + f27[dir::DIR_PMM] + f27[dir::DIR_MPM]) -
    rho - (vx * vx + vy * vy + vz * vz) * (basics::constant::c1o1 + rho)) * 
    basics::constant::c1o2 + rho; // times zero for incompressible case                 
                          // Attention: op defined directly to op = 1 ; ^^^^(1.0/op-0.5)=0.5
}


}

#endif
