//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//=======================================================================================
#ifndef LBM_INTERPOLATION_CF_H
#define LBM_INTERPOLATION_CF_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/constants/NumericConstants.h>

#include "lbm/constants/D3Q27.h"

#include "lbm/KernelParameter.h"
#include "lbm/Chimera.h"

#include "lbm/refinement/Coefficients.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

namespace vf::lbm
{

inline __host__ __device__ void interpolate_cf(real* const f, const real& omegaF, const real& eps_new, const Coefficients &coefficients, const real& x, const real& y, const real& z)
{
    const real useNEQ = c1o1;

    const real kxyAverage    = c0o1;
    const real kyzAverage    = c0o1;
    const real kxzAverage    = c0o1;
    const real kxxMyyAverage = c0o1;
    const real kxxMzzAverage = c0o1;

    const real& a_000 = coefficients.a_000;
    const real& b_000 = coefficients.b_000;
    const real& c_000 = coefficients.c_000;
    const real& d_000 = coefficients.d_000;

    const real& a_100 = coefficients.a_100;
    const real& b_100 = coefficients.b_100;
    const real& c_100 = coefficients.c_100;
    const real& d_100 = coefficients.d_100;

    const real& a_010 = coefficients.a_010;
    const real& b_010 = coefficients.b_010;
    const real& c_010 = coefficients.c_010;
    const real& d_010 = coefficients.d_010;

    const real& a_001 = coefficients.a_001;
    const real& b_001 = coefficients.b_001;
    const real& c_001 = coefficients.c_001;
    const real& d_001 = coefficients.d_001;

    const real& d_110 = coefficients.d_110, &d_101 = coefficients.d_101, &d_011 = coefficients.d_011;
    
    const real& a_200 = coefficients.a_200, &a_020 = coefficients.a_020, &a_002 = coefficients.a_002;
    const real& b_200 = coefficients.b_200, &b_020 = coefficients.b_020, &b_002 = coefficients.b_002;
    const real& c_200 = coefficients.c_200, &c_020 = coefficients.c_020, &c_002 = coefficients.c_002;

    const real& a_110 = coefficients.a_110, &a_101 = coefficients.a_101, &a_011 = coefficients.a_011;
    const real& b_110 = coefficients.b_110, &b_101 = coefficients.b_101, &b_011 = coefficients.b_011;
    const real& c_110 = coefficients.c_110, &c_101 = coefficients.c_101, &c_011 = coefficients.c_011;

    const real &a_111 = coefficients.a_111, &b_111 = coefficients.b_111, &c_111 = coefficients.c_111, &d_111 = coefficients.d_111;

    const real &LaplaceRho = coefficients.LaplaceRho;


    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set all moments to zero
    //!      
    real m_111 = c0o1;
    real m_211 = c0o1;
    real m_011 = c0o1;
    real m_121 = c0o1;
    real m_101 = c0o1;
    real m_112 = c0o1;
    real m_110 = c0o1;
    real m_221 = c0o1;
    real m_001 = c0o1;
    real m_201 = c0o1;
    real m_021 = c0o1;
    real m_212 = c0o1;
    real m_010 = c0o1;
    real m_210 = c0o1;
    real m_012 = c0o1;
    real m_122 = c0o1;
    real m_100 = c0o1;
    real m_120 = c0o1;
    real m_102 = c0o1;
    real m_222 = c0o1;
    real m_022 = c0o1;
    real m_202 = c0o1;
    real m_002 = c0o1;
    real m_220 = c0o1;
    real m_020 = c0o1;
    real m_200 = c0o1;
    real m_000 = c0o1;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Define aliases to use the same variable for the distributions (f's):
    //!
    real& f_000 = m_111;
    real& f_P00 = m_211;
    real& f_M00 = m_011;
    real& f_0P0 = m_121;
    real& f_0M0 = m_101;
    real& f_00P = m_112;
    real& f_00M = m_110;
    real& f_PP0 = m_221;
    real& f_MM0 = m_001;
    real& f_PM0 = m_201;
    real& f_MP0 = m_021;
    real& f_P0P = m_212;
    real& f_M0M = m_010;
    real& f_P0M = m_210;
    real& f_M0P = m_012;
    real& f_0PP = m_122;
    real& f_0MM = m_100;
    real& f_0PM = m_120;
    real& f_0MP = m_102;
    real& f_PPP = m_222;
    real& f_MPP = m_022;
    real& f_PMP = m_202;
    real& f_MMP = m_002;
    real& f_PPM = m_220;
    real& f_MPM = m_020;
    real& f_PMM = m_200;
    real& f_MMM = m_000;


    ////////////////////////////////////////////////////////////////////////////////
    //! - Set macroscopic values on destination node (zeroth and first order moments)
    //!
    real press = d_000 + x * d_100 + y * d_010 + z * d_001 +
                 x * y * d_110 + x * z * d_101 + y * z * d_011 + x * y * z * d_111 + c3o1 * x * x * LaplaceRho;
    real vvx   = a_000 + x * a_100 + y * a_010 + z * a_001 +
                 x * x * a_200 + y * y * a_020 + z * z * a_002 +
                 x * y * a_110 + x * z * a_101 + y * z * a_011 + x * y * z * a_111;
    real vvy   = b_000 + x * b_100 + y * b_010 + z * b_001 +
                 x * x * b_200 + y * y * b_020 + z * z * b_002 +
                 x * y * b_110 + x * z * b_101 + y * z * b_011 + x * y * z * b_111;
    real vvz   = c_000 + x * c_100 + y * c_010 + z * c_001 +
                 x * x * c_200 + y * y * c_020 + z * z * c_002 +
                 x * y * c_110 + x * z * c_101 + y * z * c_011 + x * y * z * c_111;

    m_000 = press; // m_000 is press, if drho is interpolated directly

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (second to sixth order) on destination node
    //!
    // linear combinations for second order moments
    real mxxPyyPzz = m_000;

    real mxxMyy = -c2o3 * (a_100 - b_010 + kxxMyyAverage + c2o1 * a_200 * x - b_110 * x + a_110 * y
                  -c2o1 * b_020 * y + a_101 * z - b_011 * z - b_111 * x * z + a_111 * y * z) * eps_new/ omegaF * (c1o1 + press);
    real mxxMzz = -c2o3 * (a_100 - c_001 + kxxMzzAverage + c2o1 * a_200 * x - c_101 * x + a_110 * y
                  -c_011 * y - c_111 * x * y + a_101 * z - c2o1 * c_002 * z + a_111 * y * z) * eps_new/ omegaF * (c1o1 + press);

    m_011 = -c1o3 * (b_001 + c_010 + kyzAverage + b_101 * x + c_110 * x + b_011 * y + c2o1 * c_020 * y
            + b_111 * x * y + c2o1 * b_002 * z + c_011 * z + c_111 * x * z) * eps_new / omegaF * (c1o1 + press);
    m_101 = -c1o3 * (a_001 + c_100 + kxzAverage + a_101 * x + c2o1 * c_200 * x + a_011 * y + c_110 * y
            + a_111 * x * y + c2o1 * a_002 * z + c_101 * z + c_111 * y * z) * eps_new / omegaF * (c1o1 + press);
    m_110 = -c1o3 * (a_010 + b_100 + kxyAverage + a_110 * x + c2o1 * b_200 * x + c2o1 * a_020 * y
            + b_110 * y + a_011 * z + b_101 * z + a_111 * x * z + b_111 * y * z) * eps_new / omegaF * (c1o1 + press);

    m_200 = c1o3 * (        mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m_020 = c1o3 * (-c2o1 * mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m_002 = c1o3 * (        mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * useNEQ;

    // linear combinations for third order moments
    m_111 = c0o1;

    real mxxyPyzz = c0o1;
    real mxxyMyzz = c0o1;
    real mxxzPyyz = c0o1;
    real mxxzMyyz = c0o1;
    real mxyyPxzz = c0o1;
    real mxyyMxzz = c0o1;

    m_210 = ( mxxyMyzz + mxxyPyzz) * c1o2;
    m_012 = (-mxxyMyzz + mxxyPyzz) * c1o2;
    m_201 = ( mxxzMyyz + mxxzPyyz) * c1o2;
    m_021 = (-mxxzMyyz + mxxzPyyz) * c1o2;
    m_120 = ( mxyyMxzz + mxyyPxzz) * c1o2;
    m_102 = (-mxyyMxzz + mxyyPxzz) * c1o2;

    // fourth order moments
    m_022 = m_000 * c1o9;
    m_202 = m_022;
    m_220 = m_022;

    // fifth order moments

    // sixth order moment
    m_222 = m_000 * c1o27;

    real vx_sq = vvx * vvx;
    real vy_sq = vvy * vvy;
    real vz_sq = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!

    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardInverseChimeraWithK(m_000, m_100, m_200, vvx, vx_sq, c1o1, c1o1);
    backwardChimera(            m_010, m_110, m_210, vvx, vx_sq);
    backwardInverseChimeraWithK(m_020, m_120, m_220, vvx, vx_sq, c3o1, c1o3);
    backwardChimera(            m_001, m_101, m_201, vvx, vx_sq);
    backwardChimera(            m_011, m_111, m_211, vvx, vx_sq);
    backwardChimera(            m_021, m_121, m_221, vvx, vx_sq);
    backwardInverseChimeraWithK(m_002, m_102, m_202, vvx, vx_sq, c3o1, c1o3);
    backwardChimera(            m_012, m_112, m_212, vvx, vx_sq);
    backwardInverseChimeraWithK(m_022, m_122, m_222, vvx, vx_sq, c9o1, c1o9);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardInverseChimeraWithK(m_000, m_010, m_020, vvy, vy_sq, c6o1, c1o6);
    backwardChimera(            m_001, m_011, m_021, vvy, vy_sq);
    backwardInverseChimeraWithK(m_002, m_012, m_022, vvy, vy_sq, c18o1, c1o18);
    backwardInverseChimeraWithK(m_100, m_110, m_120, vvy, vy_sq, c3o2, c2o3);
    backwardChimera(            m_101, m_111, m_121, vvy, vy_sq);
    backwardInverseChimeraWithK(m_102, m_112, m_122, vvy, vy_sq, c9o2, c2o9);
    backwardInverseChimeraWithK(m_200, m_210, m_220, vvy, vy_sq, c6o1, c1o6);
    backwardChimera(            m_201, m_211, m_221, vvy, vy_sq);
    backwardInverseChimeraWithK(m_202, m_212, m_222, vvy, vy_sq, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardInverseChimeraWithK(m_000, m_001, m_002, vvz, vz_sq, c36o1, c1o36);
    backwardInverseChimeraWithK(m_010, m_011, m_012, vvz, vz_sq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m_020, m_021, m_022, vvz, vz_sq, c36o1, c1o36);
    backwardInverseChimeraWithK(m_100, m_101, m_102, vvz, vz_sq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m_110, m_111, m_112, vvz, vz_sq, c9o4,  c4o9);
    backwardInverseChimeraWithK(m_120, m_121, m_122, vvz, vz_sq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m_200, m_201, m_202, vvz, vz_sq, c36o1, c1o36);
    backwardInverseChimeraWithK(m_210, m_211, m_212, vvz, vz_sq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m_220, m_221, m_222, vvz, vz_sq, c36o1, c1o36);

    f[DIR_000] = f_000;
    f[DIR_P00] = f_P00;
    f[DIR_M00] = f_M00;
    f[DIR_0P0] = f_0P0;
    f[DIR_0M0] = f_0M0;
    f[DIR_00P] = f_00P;
    f[DIR_00M] = f_00M;
    f[DIR_PP0] = f_PP0;
    f[DIR_MM0] = f_MM0;
    f[DIR_PM0] = f_PM0;
    f[DIR_MP0] = f_MP0;
    f[DIR_P0P] = f_P0P;
    f[DIR_M0M] = f_M0M;
    f[DIR_P0M] = f_P0M;
    f[DIR_M0P] = f_M0P;
    f[DIR_0PP] = f_0PP;
    f[DIR_0MM] = f_0MM;
    f[DIR_0PM] = f_0PM;
    f[DIR_0MP] = f_0MP;
    f[DIR_PPP] = f_PPP;
    f[DIR_MPP] = f_MPP;
    f[DIR_PMP] = f_PMP;
    f[DIR_MMP] = f_MMP;
    f[DIR_PPM] = f_PPM;
    f[DIR_MPM] = f_MPM;
    f[DIR_PMM] = f_PMM;
    f[DIR_MMM] = f_MMM;
}



}

#endif
