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
//! \addtogroup collision
//! \ingroup lbm
//! \{

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include <cuda_helper/CudaIndexCalculation.h>

#include <lbm/ChimeraTransformation.h>
#include <lbm/MacroscopicQuantities.h>
#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/constants/D3Q27.h>

#include "CollisionParameter.h"


namespace vf::lbm::advection_diffusion
{

constexpr void runF16AdvectionDiffusion(ADCollisionParameter& parameters)
{
    using namespace vf::basics::constant;
    using namespace dir;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set local distributions Advection Diffusion
    //!
    real f000 = parameters.distribution[d000];
    real fP00 = parameters.distribution[dP00];
    real fM00 = parameters.distribution[dM00];
    real f0P0 = parameters.distribution[d0P0];
    real f0M0 = parameters.distribution[d0M0];
    real f00P = parameters.distribution[d00P];
    real f00M = parameters.distribution[d00M];
    real fPP0 = parameters.distribution[dPP0];
    real fMM0 = parameters.distribution[dMM0];
    real fPM0 = parameters.distribution[dPM0];
    real fMP0 = parameters.distribution[dMP0];
    real fP0P = parameters.distribution[dP0P];
    real fM0M = parameters.distribution[dM0M];
    real fP0M = parameters.distribution[dP0M];
    real fM0P = parameters.distribution[dM0P];
    real f0PP = parameters.distribution[d0PP];
    real f0MM = parameters.distribution[d0MM];
    real f0PM = parameters.distribution[d0PM];
    real f0MP = parameters.distribution[d0MP];
    real fPPP = parameters.distribution[dPPP];
    real fMPP = parameters.distribution[dMPP];
    real fPMP = parameters.distribution[dPMP];
    real fMMP = parameters.distribution[dMMP];
    real fPPM = parameters.distribution[dPPM];
    real fMPM = parameters.distribution[dMPM];
    real fPMM = parameters.distribution[dPMM];
    real fMMM = parameters.distribution[dMMM];

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Define aliases to use the same variable for the moments (m's):
    //!
    real& m211 = fP00;
    real& m011 = fM00;
    real& m121 = f0P0;
    real& m101 = f0M0;
    real& m112 = f00P;
    real& m110 = f00M;
    real& m221 = fPP0;
    real& m001 = fMM0;
    real& m201 = fPM0;
    real& m021 = fMP0;
    real& m212 = fP0P;
    real& m010 = fM0M;
    real& m210 = fP0M;
    real& m012 = fM0P;
    real& m122 = f0PP;
    real& m100 = f0MM;
    real& m120 = f0PM;
    real& m102 = f0MP;
    real& m111 = f000;
    real& m222 = fPPP;
    real& m002 = fMMP;
    real& m202 = fPMP;
    real& m022 = fMPP;
    real& m220 = fPPM;
    real& m000 = fMMM;
    real& m200 = fPMM;
    real& m020 = fMPM;

    const real concentration = vf::lbm::getDensity(parameters.distribution);

    ////////////////////////////////////////////////////////////////////////////////////
    // calculate the square of velocities for this lattice node
    const real vvx = parameters.velocityX;
    const real vvy = parameters.velocityY;
    const real vvz = parameters.velocityZ;

    const real vx2 = vvx * vvx;
    const real vy2 = vvy * vvy;
    const real vz2 = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from distributions to central moments as defined in Eq. (43)-(45) in \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    forwardChimera(m000, m001, m002, vvz, vz2);
    forwardChimera(m010, m011, m012, vvz, vz2);
    forwardChimera(m020, m021, m022, vvz, vz2);
    forwardChimera(m100, m101, m102, vvz, vz2);
    forwardChimera(m110, m111, m112, vvz, vz2);
    forwardChimera(m120, m121, m122, vvz, vz2);
    forwardChimera(m200, m201, m202, vvz, vz2);
    forwardChimera(m210, m211, m212, vvz, vz2);
    forwardChimera(m220, m221, m222, vvz, vz2);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    forwardChimera(m000, m010, m020, vvy, vy2);
    forwardChimera(m001, m011, m021, vvy, vy2);
    forwardChimera(m002, m012, m022, vvy, vy2);
    forwardChimera(m100, m110, m120, vvy, vy2);
    forwardChimera(m101, m111, m121, vvy, vy2);
    forwardChimera(m102, m112, m122, vvy, vy2);
    forwardChimera(m200, m210, m220, vvy, vy2);
    forwardChimera(m201, m211, m221, vvy, vy2);
    forwardChimera(m202, m212, m222, vvy, vy2);

    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    forwardChimera(m000, m100, m200, vvx, vx2);
    forwardChimera(m010, m110, m210, vvx, vx2);
    forwardChimera(m020, m120, m220, vvx, vx2);
    forwardChimera(m001, m101, m201, vvx, vx2);
    forwardChimera(m011, m111, m211, vvx, vx2);
    forwardChimera(m021, m121, m221, vvx, vx2);
    forwardChimera(m002, m102, m202, vvx, vx2);
    forwardChimera(m012, m112, m212, vvx, vx2);
    forwardChimera(m022, m122, m222, vvx, vx2);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Factorized central moments for Advection Diffusion Equation - Eq. (15)-(16) in \ref
    //! <a href="https://doi.org/10.1016/j.advwatres.2015.09.015"><b>[ X. Yang et al. (2016),
    //! DOI: 10.1016/j.advwatres.2015.09.015]</b></a>
    //!

    // linearized orthogonalization of 3rd order central moments
    real M012 = m012 - m010 * c1o3;
    real M120 = m120 - m100 * c1o3;
    real M021 = m021 - m001 * c1o3;
    real M210 = m210 - m010 * c1o3;
    real M201 = m201 - m001 * c1o3;
    real M102 = m102 - m100 * c1o3;
    // linearized orthogonalization of 5th order central moments
    real M212 = m212 - m010 * c1o9;
    real M122 = m122 - m100 * c1o9;
    real M221 = m221 - m001 * c1o9;

    // collision of 1st order moments
    m100 *= c1o1 - parameters.omega;
    m010 *= c1o1 - parameters.omega;
    m001 *= c1o1 - parameters.omega;

    // equilibration of 3rd order moments
    M012 = c0o1;
    M120 = c0o1;
    M021 = c0o1;
    M210 = c0o1;
    M201 = c0o1;
    M102 = c0o1;
    m111 = c0o1;

    // equilibration of 5th order moments
    M212 = c0o1;
    M122 = c0o1;
    M221 = c0o1;

    // equilibration of 2nd order moments
    m110 = c0o1;
    m101 = c0o1;
    m011 = c0o1;

    m200 = c1o3 * concentration;
    m020 = c1o3 * concentration;
    m002 = c1o3 * concentration;

    // equilibration of 4th order moments
    m022 = c1o9 * concentration;
    m202 = c1o9 * concentration;
    m220 = c1o9 * concentration;

    m211 = c0o1;
    m121 = c0o1;
    m112 = c0o1;

    // equilibration of 6th order moment
    m222 = c1o27 * concentration;

    // from linearized orthogonalization 3rd order central moments to central moments
    m012 = M012 + m010 * c1o3;
    m120 = M120 + m100 * c1o3;
    m021 = M021 + m001 * c1o3;
    m210 = M210 + m010 * c1o3;
    m201 = M201 + m001 * c1o3;
    m102 = M102 + m100 * c1o3;

    // from linearized orthogonalization 5th order central moments to central moments
    m212 = M212 + m010 * c1o9;
    m122 = M122 + m100 * c1o9;
    m221 = M221 + m001 * c1o9;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from  central moments to distributions as defined in Eq. (88)-(96) in \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardChimera(m000, m100, m200, vvx, vx2);
    backwardChimera(m010, m110, m210, vvx, vx2);
    backwardChimera(m020, m120, m220, vvx, vx2);
    backwardChimera(m001, m101, m201, vvx, vx2);
    backwardChimera(m011, m111, m211, vvx, vx2);
    backwardChimera(m021, m121, m221, vvx, vx2);
    backwardChimera(m002, m102, m202, vvx, vx2);
    backwardChimera(m012, m112, m212, vvx, vx2);
    backwardChimera(m022, m122, m222, vvx, vx2);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardChimera(m000, m010, m020, vvy, vy2);
    backwardChimera(m001, m011, m021, vvy, vy2);
    backwardChimera(m002, m012, m022, vvy, vy2);
    backwardChimera(m100, m110, m120, vvy, vy2);
    backwardChimera(m101, m111, m121, vvy, vy2);
    backwardChimera(m102, m112, m122, vvy, vy2);
    backwardChimera(m200, m210, m220, vvy, vy2);
    backwardChimera(m201, m211, m221, vvy, vy2);
    backwardChimera(m202, m212, m222, vvy, vy2);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardChimera(m000, m001, m002, vvz, vz2);
    backwardChimera(m010, m011, m012, vvz, vz2);
    backwardChimera(m020, m021, m022, vvz, vz2);
    backwardChimera(m100, m101, m102, vvz, vz2);
    backwardChimera(m110, m111, m112, vvz, vz2);
    backwardChimera(m120, m121, m122, vvz, vz2);
    backwardChimera(m200, m201, m202, vvz, vz2);
    backwardChimera(m210, m211, m212, vvz, vz2);
    backwardChimera(m220, m221, m222, vvz, vz2);

    parameters.distribution[dP00] = fP00;
    parameters.distribution[dM00] = fM00;
    parameters.distribution[d0P0] = f0P0;
    parameters.distribution[d0M0] = f0M0;
    parameters.distribution[d00P] = f00P;
    parameters.distribution[d00M] = f00M;
    parameters.distribution[dPP0] = fPP0;
    parameters.distribution[dMM0] = fMM0;
    parameters.distribution[dPM0] = fPM0;
    parameters.distribution[dMP0] = fMP0;
    parameters.distribution[dP0P] = fP0P;
    parameters.distribution[dM0M] = fM0M;
    parameters.distribution[dP0M] = fP0M;
    parameters.distribution[dM0P] = fM0P;
    parameters.distribution[d0PP] = f0PP;
    parameters.distribution[d0MM] = f0MM;
    parameters.distribution[d0PM] = f0PM;
    parameters.distribution[d0MP] = f0MP;
    parameters.distribution[d000] = f000;
    parameters.distribution[dPPP] = fPPP;
    parameters.distribution[dPMP] = fPMP;
    parameters.distribution[dPPM] = fPPM;
    parameters.distribution[dPMM] = fPMM;
    parameters.distribution[dMPP] = fMPP;
    parameters.distribution[dMMP] = fMMP;
    parameters.distribution[dMPM] = fMPM;
    parameters.distribution[dMMM] = fMMM;

    parameters.concentration = concentration;
}
////////////////////////////////////////////////////////////////////////////////

}
//! \}
