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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup refinement
//! \ingroup lbm
//! \{
//=======================================================================================
#ifndef LBM_INTERPOLATION_FC_H
#define LBM_INTERPOLATION_FC_H

#include <basics/constants/NumericConstants.h>

#include "lbm/constants/D3Q27.h"
#include "lbm/ChimeraTransformation.h"

#include "lbm/interpolation/InterpolationCoefficients.h"

namespace vf::lbm
{

constexpr void interpolateFC(real* const f, const real epsnew, const real omegaC, const InterpolationCoefficients& coefficients)
{
    using namespace vf::basics::constant;

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

    constexpr real useNEQ = c1o1;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set macroscopic values on destination node (zeroth and first order moments)
    //!
    const real press = coefficients.d000 - c2o1 * coefficients.LaplaceRho * c1o8;
    const real vvx = coefficients.a000;
    const real vvy = coefficients.b000;
    const real vvz = coefficients.c000;

    m000 = press; // m000 is press, if drho is interpolated directly

    const real vxsq = vvx * vvx;
    const real vysq = vvy * vvy;
    const real vzsq = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (second to sixth order) on destination node
    //!
    // linear combinations for second order moments
    const real mxxPyyPzz = m000;

    constexpr real kxyAverage = c0o1;
    constexpr real kyzAverage = c0o1;
    constexpr real kxzAverage = c0o1;
    constexpr real kxxMyyAverage = c0o1;
    constexpr real kxxMzzAverage = c0o1;

    const real mxxMyy = -c2o3 * ((coefficients.a100 - coefficients.b010) + kxxMyyAverage) * epsnew / omegaC * (c1o1 + press);
    const real mxxMzz = -c2o3 * ((coefficients.a100 - coefficients.c001) + kxxMzzAverage) * epsnew / omegaC * (c1o1 + press);

    m011 = -c1o3 * ((coefficients.b001 + coefficients.c010) + kyzAverage) * epsnew / omegaC * (c1o1 + press);
    m101 = -c1o3 * ((coefficients.a001 + coefficients.c100) + kxzAverage) * epsnew / omegaC * (c1o1 + press);
    m110 = -c1o3 * ((coefficients.a010 + coefficients.b100) + kxyAverage) * epsnew / omegaC * (c1o1 + press);

    m200 = c1o3 * (        mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m020 = c1o3 * (-c2o1 * mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m002 = c1o3 * (        mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * useNEQ;

    // linear combinations for third order moments
    m111 = c0o1;

    constexpr real mxxyPyzz = c0o1;
    constexpr real mxxyMyzz = c0o1;
    constexpr real mxxzPyyz = c0o1;
    constexpr real mxxzMyyz = c0o1;
    constexpr real mxyyPxzz = c0o1;
    constexpr real mxyyMxzz = c0o1;

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

    // sixth order moments
    m222 = m000 * c1o27;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardChimeraWithInverseK(m000, m100, m200, vvx, vxsq, c1o1, c1o1);
    backwardChimera(            m010, m110, m210, vvx, vxsq);
    backwardChimeraWithInverseK(m020, m120, m220, vvx, vxsq, c3o1, c1o3);
    backwardChimera(            m001, m101, m201, vvx, vxsq);
    backwardChimera(            m011, m111, m211, vvx, vxsq);
    backwardChimera(            m021, m121, m221, vvx, vxsq);
    backwardChimeraWithInverseK(m002, m102, m202, vvx, vxsq, c3o1, c1o3);
    backwardChimera(            m012, m112, m212, vvx, vxsq);
    backwardChimeraWithInverseK(m022, m122, m222, vvx, vxsq, c9o1, c1o9);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardChimeraWithInverseK(m000, m010, m020, vvy, vysq, c6o1, c1o6);
    backwardChimera(            m001, m011, m021, vvy, vysq);
    backwardChimeraWithInverseK(m002, m012, m022, vvy, vysq, c18o1, c1o18);
    backwardChimeraWithInverseK(m100, m110, m120, vvy, vysq, c3o2, c2o3);
    backwardChimera(            m101, m111, m121, vvy, vysq);
    backwardChimeraWithInverseK(m102, m112, m122, vvy, vysq, c9o2, c2o9);
    backwardChimeraWithInverseK(m200, m210, m220, vvy, vysq, c6o1, c1o6);
    backwardChimera(            m201, m211, m221, vvy, vysq);
    backwardChimeraWithInverseK(m202, m212, m222, vvy, vysq, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardChimeraWithInverseK(m000, m001, m002, vvz, vzsq, c36o1, c1o36);
    backwardChimeraWithInverseK(m010, m011, m012, vvz, vzsq, c9o1,  c1o9);
    backwardChimeraWithInverseK(m020, m021, m022, vvz, vzsq, c36o1, c1o36);
    backwardChimeraWithInverseK(m100, m101, m102, vvz, vzsq, c9o1,  c1o9);
    backwardChimeraWithInverseK(m110, m111, m112, vvz, vzsq, c9o4,  c4o9);
    backwardChimeraWithInverseK(m120, m121, m122, vvz, vzsq, c9o1,  c1o9);
    backwardChimeraWithInverseK(m200, m201, m202, vvz, vzsq, c36o1, c1o36);
    backwardChimeraWithInverseK(m210, m211, m212, vvz, vzsq, c9o1,  c1o9);
    backwardChimeraWithInverseK(m220, m221, m222, vvz, vzsq, c36o1, c1o36);

    f[dir::d000] = f000;
    f[dir::dP00] = fP00;
    f[dir::dM00] = fM00;
    f[dir::d0P0] = f0P0;
    f[dir::d0M0] = f0M0;
    f[dir::d00P] = f00P;
    f[dir::d00M] = f00M;
    f[dir::dPP0] = fPP0;
    f[dir::dMM0] = fMM0;
    f[dir::dPM0] = fPM0;
    f[dir::dMP0] = fMP0;
    f[dir::dP0P] = fP0P;
    f[dir::dM0M] = fM0M;
    f[dir::dP0M] = fP0M;
    f[dir::dM0P] = fM0P;
    f[dir::d0PP] = f0PP;
    f[dir::d0MM] = f0MM;
    f[dir::d0PM] = f0PM;
    f[dir::d0MP] = f0MP;
    f[dir::dPPP] = fPPP;
    f[dir::dMPP] = fMPP;
    f[dir::dPMP] = fPMP;
    f[dir::dMMP] = fMMP;
    f[dir::dPPM] = fPPM;
    f[dir::dMPM] = fMPM;
    f[dir::dPMM] = fPMM;
    f[dir::dMMM] = fMMM;
}

}

#endif

//! \}
