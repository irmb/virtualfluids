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

#include "lbm/ChimeraTransformation.h"

#include "lbm/interpolation/InterpolationCoefficients.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

namespace vf::lbm
{

inline __host__ __device__ void interpolateCF(real* const f, const real& omegaF, const real& epsnew, const InterpolationCoefficients &coefficients, const real& x, const real& y, const real& z)
{
    const real useNEQ = c1o1;

    const real kxyAverage    = c0o1;
    const real kyzAverage    = c0o1;
    const real kxzAverage    = c0o1;
    const real kxxMyyAverage = c0o1;
    const real kxxMzzAverage = c0o1;

    const real& a000 = coefficients.a000;
    const real& b000 = coefficients.b000;
    const real& c000 = coefficients.c000;
    const real& d000 = coefficients.d000;

    const real& a100 = coefficients.a100;
    const real& b100 = coefficients.b100;
    const real& c100 = coefficients.c100;
    const real& d100 = coefficients.d100;

    const real& a010 = coefficients.a010;
    const real& b010 = coefficients.b010;
    const real& c010 = coefficients.c010;
    const real& d010 = coefficients.d010;

    const real& a001 = coefficients.a001;
    const real& b001 = coefficients.b001;
    const real& c001 = coefficients.c001;
    const real& d001 = coefficients.d001;

    const real& d110 = coefficients.d110, &d101 = coefficients.d101, &d011 = coefficients.d011;
    
    const real& a200 = coefficients.a200, &a020 = coefficients.a020, &a002 = coefficients.a002;
    const real& b200 = coefficients.b200, &b020 = coefficients.b020, &b002 = coefficients.b002;
    const real& c200 = coefficients.c200, &c020 = coefficients.c020, &c002 = coefficients.c002;

    const real& a110 = coefficients.a110, &a101 = coefficients.a101, &a011 = coefficients.a011;
    const real& b110 = coefficients.b110, &b101 = coefficients.b101, &b011 = coefficients.b011;
    const real& c110 = coefficients.c110, &c101 = coefficients.c101, &c011 = coefficients.c011;

    const real &a111 = coefficients.a111, &b111 = coefficients.b111, &c111 = coefficients.c111, &d111 = coefficients.d111;

    const real &LaplaceRho = coefficients.LaplaceRho;


    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set all moments to zero
    //!      
    real m111 = c0o1;
    real m211 = c0o1;
    real m011 = c0o1;
    real m121 = c0o1;
    real m101 = c0o1;
    real m112 = c0o1;
    real m110 = c0o1;
    real m221 = c0o1;
    real m001 = c0o1;
    real m201 = c0o1;
    real m021 = c0o1;
    real m212 = c0o1;
    real m010 = c0o1;
    real m210 = c0o1;
    real m012 = c0o1;
    real m122 = c0o1;
    real m100 = c0o1;
    real m120 = c0o1;
    real m102 = c0o1;
    real m222 = c0o1;
    real m022 = c0o1;
    real m202 = c0o1;
    real m002 = c0o1;
    real m220 = c0o1;
    real m020 = c0o1;
    real m200 = c0o1;
    real m000 = c0o1;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Define aliases to use the same variable for the distributions (f's):
    //!
    real& f000 = m111;
    real& fP00 = m211;
    real& fM00 = m011;
    real& f0P0 = m121;
    real& f0M0 = m101;
    real& f00P = m112;
    real& f00M = m110;
    real& fPP0 = m221;
    real& fMM0 = m001;
    real& fPM0 = m201;
    real& fMP0 = m021;
    real& fP0P = m212;
    real& fM0M = m010;
    real& fP0M = m210;
    real& fM0P = m012;
    real& f0PP = m122;
    real& f0MM = m100;
    real& f0PM = m120;
    real& f0MP = m102;
    real& fPPP = m222;
    real& fMPP = m022;
    real& fPMP = m202;
    real& fMMP = m002;
    real& fPPM = m220;
    real& fMPM = m020;
    real& fPMM = m200;
    real& fMMM = m000;


    ////////////////////////////////////////////////////////////////////////////////
    //! - Set macroscopic values on destination node (zeroth and first order moments)
    //!
    real press = d000 + x * d100 + y * d010 + z * d001 +
                 x * y * d110 + x * z * d101 + y * z * d011 + x * y * z * d111 + c3o1 * x * x * LaplaceRho;
    real vvx   = a000 + x * a100 + y * a010 + z * a001 +
                 x * x * a200 + y * y * a020 + z * z * a002 +
                 x * y * a110 + x * z * a101 + y * z * a011 + x * y * z * a111;
    real vvy   = b000 + x * b100 + y * b010 + z * b001 +
                 x * x * b200 + y * y * b020 + z * z * b002 +
                 x * y * b110 + x * z * b101 + y * z * b011 + x * y * z * b111;
    real vvz   = c000 + x * c100 + y * c010 + z * c001 +
                 x * x * c200 + y * y * c020 + z * z * c002 +
                 x * y * c110 + x * z * c101 + y * z * c011 + x * y * z * c111;

    m000 = press; // m000 is press, if drho is interpolated directly

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (second to sixth order) on destination node
    //!
    // linear combinations for second order moments
    real mxxPyyPzz = m000;

    real mxxMyy = -c2o3 * (a100 - b010 + kxxMyyAverage + c2o1 * a200 * x - b110 * x + a110 * y
                  -c2o1 * b020 * y + a101 * z - b011 * z - b111 * x * z + a111 * y * z) * epsnew/ omegaF * (c1o1 + press);
    real mxxMzz = -c2o3 * (a100 - c001 + kxxMzzAverage + c2o1 * a200 * x - c101 * x + a110 * y
                  -c011 * y - c111 * x * y + a101 * z - c2o1 * c002 * z + a111 * y * z) * epsnew/ omegaF * (c1o1 + press);

    m011 = -c1o3 * (b001 + c010 + kyzAverage + b101 * x + c110 * x + b011 * y + c2o1 * c020 * y
            + b111 * x * y + c2o1 * b002 * z + c011 * z + c111 * x * z) * epsnew / omegaF * (c1o1 + press);
    m101 = -c1o3 * (a001 + c100 + kxzAverage + a101 * x + c2o1 * c200 * x + a011 * y + c110 * y
            + a111 * x * y + c2o1 * a002 * z + c101 * z + c111 * y * z) * epsnew / omegaF * (c1o1 + press);
    m110 = -c1o3 * (a010 + b100 + kxyAverage + a110 * x + c2o1 * b200 * x + c2o1 * a020 * y
            + b110 * y + a011 * z + b101 * z + a111 * x * z + b111 * y * z) * epsnew / omegaF * (c1o1 + press);

    m200 = c1o3 * (        mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m020 = c1o3 * (-c2o1 * mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m002 = c1o3 * (        mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * useNEQ;

    // linear combinations for third order moments
    m111 = c0o1;

    real mxxyPyzz = c0o1;
    real mxxyMyzz = c0o1;
    real mxxzPyyz = c0o1;
    real mxxzMyyz = c0o1;
    real mxyyPxzz = c0o1;
    real mxyyMxzz = c0o1;

    m210 = ( mxxyMyzz + mxxyPyzz) * c1o2;
    m012 = (-mxxyMyzz + mxxyPyzz) * c1o2;
    m201 = ( mxxzMyyz + mxxzPyyz) * c1o2;
    m021 = (-mxxzMyyz + mxxzPyyz) * c1o2;
    m120 = ( mxyyMxzz + mxyyPxzz) * c1o2;
    m102 = (-mxyyMxzz + mxyyPxzz) * c1o2;

    // fourth order moments
    m022 = m000 * c1o9;
    m202 = m022;
    m220 = m022;

    // fifth order moments

    // sixth order moment
    m222 = m000 * c1o27;

    real vxsq = vvx * vvx;
    real vysq = vvy * vvy;
    real vzsq = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!

    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardInverseChimeraWithK(m000, m100, m200, vvx, vxsq, c1o1, c1o1);
    backwardChimera(            m010, m110, m210, vvx, vxsq);
    backwardInverseChimeraWithK(m020, m120, m220, vvx, vxsq, c3o1, c1o3);
    backwardChimera(            m001, m101, m201, vvx, vxsq);
    backwardChimera(            m011, m111, m211, vvx, vxsq);
    backwardChimera(            m021, m121, m221, vvx, vxsq);
    backwardInverseChimeraWithK(m002, m102, m202, vvx, vxsq, c3o1, c1o3);
    backwardChimera(            m012, m112, m212, vvx, vxsq);
    backwardInverseChimeraWithK(m022, m122, m222, vvx, vxsq, c9o1, c1o9);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardInverseChimeraWithK(m000, m010, m020, vvy, vysq, c6o1, c1o6);
    backwardChimera(            m001, m011, m021, vvy, vysq);
    backwardInverseChimeraWithK(m002, m012, m022, vvy, vysq, c18o1, c1o18);
    backwardInverseChimeraWithK(m100, m110, m120, vvy, vysq, c3o2, c2o3);
    backwardChimera(            m101, m111, m121, vvy, vysq);
    backwardInverseChimeraWithK(m102, m112, m122, vvy, vysq, c9o2, c2o9);
    backwardInverseChimeraWithK(m200, m210, m220, vvy, vysq, c6o1, c1o6);
    backwardChimera(            m201, m211, m221, vvy, vysq);
    backwardInverseChimeraWithK(m202, m212, m222, vvy, vysq, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardInverseChimeraWithK(m000, m001, m002, vvz, vzsq, c36o1, c1o36);
    backwardInverseChimeraWithK(m010, m011, m012, vvz, vzsq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m020, m021, m022, vvz, vzsq, c36o1, c1o36);
    backwardInverseChimeraWithK(m100, m101, m102, vvz, vzsq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m110, m111, m112, vvz, vzsq, c9o4,  c4o9);
    backwardInverseChimeraWithK(m120, m121, m122, vvz, vzsq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m200, m201, m202, vvz, vzsq, c36o1, c1o36);
    backwardInverseChimeraWithK(m210, m211, m212, vvz, vzsq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m220, m221, m222, vvz, vzsq, c36o1, c1o36);

    f[dir::d000] = f000;
    f[dP00] = fP00;
    f[dM00] = fM00;
    f[d0P0] = f0P0;
    f[d0M0] = f0M0;
    f[d00P] = f00P;
    f[d00M] = f00M;
    f[dPP0] = fPP0;
    f[dMM0] = fMM0;
    f[dPM0] = fPM0;
    f[dMP0] = fMP0;
    f[dP0P] = fP0P;
    f[dM0M] = fM0M;
    f[dP0M] = fP0M;
    f[dM0P] = fM0P;
    f[d0PP] = f0PP;
    f[d0MM] = f0MM;
    f[d0PM] = f0PM;
    f[d0MP] = f0MP;
    f[dPPP] = fPPP;
    f[dMPP] = fMPP;
    f[dPMP] = fPMP;
    f[dMMP] = fMMP;
    f[dPPM] = fPPM;
    f[dMPM] = fMPM;
    f[dPMM] = fPMM;
    f[dMMM] = fMMM;
}



}

#endif
