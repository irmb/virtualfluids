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
#ifndef LBM_INTERPOLATION_FC_H
#define LBM_INTERPOLATION_FC_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif


#include "lbm/constants/D3Q27.h"
#include "lbm/constants/NumericConstants.h"

#include "lbm/KernelParameter.h"
#include "lbm/Chimera.h"

#include "lbm/refinement/Coefficients.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

namespace vf::lbm
{

__host__ __device__ __inline__ void interpolate_fc(real* const f, const real eps_new, const real omegaC, Coefficients& coefficients)
{

    const real kxyAverage    = c0o1;
    const real kyzAverage    = c0o1;
    const real kxzAverage    = c0o1;
    const real kxxMyyAverage = c0o1;
    const real kxxMzzAverage = c0o1;

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
    //! - Declare local variables for destination nodes
    //!
    real vvx, vvy, vvz, vx_sq, vy_sq, vz_sq;
    real mxxPyyPzz, mxxMyy, mxxMzz, mxxyPyzz, mxxyMyzz, mxxzPyyz, mxxzMyyz, mxyyPxzz, mxyyMxzz;
    real useNEQ = c1o1; // zero; //one;   //.... one = on ..... zero = off
    real press;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set macroscopic values on destination node (zeroth and first order moments)
    //!
    press = coefficients.d_000 - c2o1 * coefficients.LaplaceRho * c1o8;
    vvx   = coefficients.a_000;
    vvy   = coefficients.b_000;
    vvz   = coefficients.c_000;

    m_000 = press; // m_000 is press, if drho is interpolated directly

    vx_sq = vvx * vvx;
    vy_sq = vvy * vvy;
    vz_sq = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (second to sixth order) on destination node
    //!
    // linear combinations for second order moments
    mxxPyyPzz = m_000;

    mxxMyy = -c2o3 * ((coefficients.a_100 - coefficients.b_010) + kxxMyyAverage) * eps_new / omegaC * (c1o1 + press);
    mxxMzz = -c2o3 * ((coefficients.a_100 - coefficients.c_001) + kxxMzzAverage) * eps_new / omegaC * (c1o1 + press);

    m_011 = -c1o3 * ((coefficients.b_001 + coefficients.c_010) + kyzAverage) * eps_new / omegaC * (c1o1 + press);
    m_101 = -c1o3 * ((coefficients.a_001 + coefficients.c_100) + kxzAverage) * eps_new / omegaC * (c1o1 + press);
    m_110 = -c1o3 * ((coefficients.a_010 + coefficients.b_100) + kxyAverage) * eps_new / omegaC * (c1o1 + press);

    m_200 = c1o3 * (        mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m_020 = c1o3 * (-c2o1 * mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m_002 = c1o3 * (        mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * useNEQ;

    // linear combinations for third order moments
    m_111 = c0o1;

    mxxyPyzz = c0o1;
    mxxyMyzz = c0o1;
    mxxzPyyz = c0o1;
    mxxzMyyz = c0o1;
    mxyyPxzz = c0o1;
    mxyyMxzz = c0o1;

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

    // sixth order moments
    m_222 = m_000 * c1o27;

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
