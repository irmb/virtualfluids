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
#ifndef ADVECTION_DIFFUSION_INTERPOLATION_FC_H
#define ADVECTION_DIFFUSION_INTERPOLATION_FC_H

#include <basics/constants/NumericConstants.h>

#include "lbm/constants/D3Q27.h"

#include "lbm/ChimeraTransformation.h"

#include "lbm/interpolation/InterpolationCoefficientsAdvectionDiffusion.h"

namespace vf::lbm::ad
{

constexpr void interpolateAdvectionDiffusionFC(real* const populations, real* populationsAD,
                                               const real omegaDiffusivityCoarse, const real epsnew,
                                               const InterpolationCoefficients& coefficients)
{
    using namespace ::vf::basics::constant;
    using namespace ::vf::lbm::dir;
    using namespace ::vf::lbm;

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
    real& populationsAD_000 = m111;
    real& populationsAD_P00 = m211;
    real& populationsAD_M00 = m011;
    real& populationsAD_0P0 = m121;
    real& populationsAD_0M0 = m101;
    real& populationsAD_00P = m112;
    real& populationsAD_00M = m110;
    real& populationsAD_PP0 = m221;
    real& populationsAD_MM0 = m001;
    real& populationsAD_PM0 = m201;
    real& populationsAD_MP0 = m021;
    real& populationsAD_P0P = m212;
    real& populationsAD_M0M = m010;
    real& populationsAD_P0M = m210;
    real& populationsAD_M0P = m012;
    real& populationsAD_0PP = m122;
    real& populationsAD_0MM = m100;
    real& populationsAD_0PM = m120;
    real& populationsAD_0MP = m102;
    real& populationsAD_PPP = m222;
    real& populationsAD_MPP = m022;
    real& populationsAD_PMP = m202;
    real& populationsAD_MMP = m002;
    real& populationsAD_PPM = m220;
    real& populationsAD_MPM = m020;
    real& populationsAD_PMM = m200;
    real& populationsAD_MMM = m000;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set macroscopic values on destination node (zeroth and first order moments)
    //!
    real press = c0o1;
    real vvx = c0o1;
    real vvy = c0o1;
    real vvz = c0o1;

    getCompressibleMacroscopicValues(populations, press, vvx, vvy, vvz);

    m000 = coefficients.d000;
    // m000 is the concentration

    real dxConcentration = coefficients.d100;
    real dyConcentration = coefficients.d010;
    real dzConcentration = coefficients.d001;

    m100 = -(c1o1) / (c3o1 * omegaDiffusivityCoarse) * c1o2 * dxConcentration;
    m010 = -(c1o1) / (c3o1 * omegaDiffusivityCoarse) * c1o2 * dyConcentration;
    m001 = -(c1o1) / (c3o1 * omegaDiffusivityCoarse) * c1o2 * dzConcentration;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (second to sixth order) on destination node
    //!
    // linear combinations for second order moments
    // const real mxxPyyPzz = m000;

    m200 = m000 * c1o3;
    m020 = m000 * c1o3;
    m002 = m000 * c1o3;

    // fourth order moments
    m022 = m000 * c1o9;
    m202 = m022;
    m220 = m022;

    // fifth order moments

    // sixth order moment
    m222 = m000 * c1o27;

    const real vxsq = vvx * vvx;
    const real vysq = vvy * vvy;
    const real vzsq = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardChimera(m000, m100, m200, vvx, vxsq);
    backwardChimera(m010, m110, m210, vvx, vxsq);
    backwardChimera(m020, m120, m220, vvx, vxsq);
    backwardChimera(m001, m101, m201, vvx, vxsq);
    backwardChimera(m011, m111, m211, vvx, vxsq);
    backwardChimera(m021, m121, m221, vvx, vxsq);
    backwardChimera(m002, m102, m202, vvx, vxsq);
    backwardChimera(m012, m112, m212, vvx, vxsq);
    backwardChimera(m022, m122, m222, vvx, vxsq);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardChimera(m000, m010, m020, vvy, vysq);
    backwardChimera(m001, m011, m021, vvy, vysq);
    backwardChimera(m002, m012, m022, vvy, vysq);
    backwardChimera(m100, m110, m120, vvy, vysq);
    backwardChimera(m101, m111, m121, vvy, vysq);
    backwardChimera(m102, m112, m122, vvy, vysq);
    backwardChimera(m200, m210, m220, vvy, vysq);
    backwardChimera(m201, m211, m221, vvy, vysq);
    backwardChimera(m202, m212, m222, vvy, vysq);
    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardChimera(m000, m001, m002, vvz, vzsq);
    backwardChimera(m010, m011, m012, vvz, vzsq);
    backwardChimera(m020, m021, m022, vvz, vzsq);
    backwardChimera(m100, m101, m102, vvz, vzsq);
    backwardChimera(m110, m111, m112, vvz, vzsq);
    backwardChimera(m120, m121, m122, vvz, vzsq);
    backwardChimera(m200, m201, m202, vvz, vzsq);
    backwardChimera(m210, m211, m212, vvz, vzsq);
    backwardChimera(m220, m221, m222, vvz, vzsq);


    populations[d000] = populationsAD_000;
    populations[dP00] = populationsAD_P00;
    populations[dM00] = populationsAD_M00;
    populations[d0P0] = populationsAD_0P0;
    populations[d0M0] = populationsAD_0M0;
    populations[d00P] = populationsAD_00P;
    populations[d00M] = populationsAD_00M;
    populations[dPP0] = populationsAD_PP0;
    populations[dMM0] = populationsAD_MM0;
    populations[dPM0] = populationsAD_PM0;
    populations[dMP0] = populationsAD_MP0;
    populations[dP0P] = populationsAD_P0P;
    populations[dM0M] = populationsAD_M0M;
    populations[dP0M] = populationsAD_P0M;
    populations[dM0P] = populationsAD_M0P;
    populations[d0PP] = populationsAD_0PP;
    populations[d0MM] = populationsAD_0MM;
    populations[d0PM] = populationsAD_0PM;
    populations[d0MP] = populationsAD_0MP;
    populations[dPPP] = populationsAD_PPP;
    populations[dMPP] = populationsAD_MPP;
    populations[dPMP] = populationsAD_PMP;
    populations[dMMP] = populationsAD_MMP;
    populations[dPPM] = populationsAD_PPM;
    populations[dMPM] = populationsAD_MPM;
    populations[dPMM] = populationsAD_PMM;
    populations[dMMM] = populationsAD_MMM;
}

} // namespace vf::lbm::ad

#endif

//! \}
