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
    auto& distribution = parameters.distributions;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set local distributions Advection Diffusion
    //!
    real& mfcbb = distribution[dP00];
    real& mfabb = distribution[dM00];
    real& mfbcb = distribution[d0P0];
    real& mfbab = distribution[d0M0];
    real& mfbbc = distribution[d00P];
    real& mfbba = distribution[d00M];
    real& mfccb = distribution[dPP0];
    real& mfaab = distribution[dMM0];
    real& mfcab = distribution[dPM0];
    real& mfacb = distribution[dMP0];
    real& mfcbc = distribution[dP0P];
    real& mfaba = distribution[dM0M];
    real& mfcba = distribution[dP0M];
    real& mfabc = distribution[dM0P];
    real& mfbcc = distribution[d0PP];
    real& mfbaa = distribution[d0MM];
    real& mfbca = distribution[d0PM];
    real& mfbac = distribution[d0MP];
    real& mfbbb = distribution[d000];
    real& mfccc = distribution[dPPP];
    real& mfaac = distribution[dMMP];
    real& mfcac = distribution[dPMP];
    real& mfacc = distribution[dMPP];
    real& mfcca = distribution[dPPM];
    real& mfaaa = distribution[dMMM];
    real& mfcaa = distribution[dPMM];
    real& mfaca = distribution[dMPM];

    const real concentration = vf::lbm::getDensity(distribution);

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
    forwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
    forwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
    forwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
    forwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
    forwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
    forwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
    forwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
    forwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
    forwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    forwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
    forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
    forwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
    forwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
    forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
    forwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
    forwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
    forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
    forwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    forwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
    forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
    forwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
    forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
    forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
    forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
    forwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
    forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
    forwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Factorized central moments for Advection Diffusion Equation - Eq. (15)-(16) in \ref
    //! <a href="https://doi.org/10.1016/j.advwatres.2015.09.015"><b>[ X. Yang et al. (2016),
    //! DOI: 10.1016/j.advwatres.2015.09.015]</b></a>
    //!

    // linearized orthogonalization of 3rd order central moments
    real Mabc = mfabc - mfaba * c1o3;
    real Mbca = mfbca - mfbaa * c1o3;
    real Macb = mfacb - mfaab * c1o3;
    real Mcba = mfcba - mfaba * c1o3;
    real Mcab = mfcab - mfaab * c1o3;
    real Mbac = mfbac - mfbaa * c1o3;
    // linearized orthogonalization of 5th order central moments
    real Mcbc = mfcbc - mfaba * c1o9;
    real Mbcc = mfbcc - mfbaa * c1o9;
    real Mccb = mfccb - mfaab * c1o9;

    // collision of 1st order moments
    mfbaa *= c1o1 - parameters.omega;
    mfaba *= c1o1 - parameters.omega;
    mfaab *= c1o1 - parameters.omega;

    // equilibration of 3rd order moments
    Mabc = c0o1;
    Mbca = c0o1;
    Macb = c0o1;
    Mcba = c0o1;
    Mcab = c0o1;
    Mbac = c0o1;
    mfbbb = c0o1;

    // equilibration of 5th order moments
    Mcbc = c0o1;
    Mbcc = c0o1;
    Mccb = c0o1;

    // equilibration of 2nd order moments
    mfbba = c0o1;
    mfbab = c0o1;
    mfabb = c0o1;

    mfcaa = c1o3 * concentration;
    mfaca = c1o3 * concentration;
    mfaac = c1o3 * concentration;

    // equilibration of 4th order moments
    mfacc = c1o9 * concentration;
    mfcac = c1o9 * concentration;
    mfcca = c1o9 * concentration;

    mfcbb = c0o1;
    mfbcb = c0o1;
    mfbbc = c0o1;

    // equilibration of 6th order moment
    mfccc = c1o27 * concentration;

    // from linearized orthogonalization 3rd order central moments to central moments
    mfabc = Mabc + mfaba * c1o3;
    mfbca = Mbca + mfbaa * c1o3;
    mfacb = Macb + mfaab * c1o3;
    mfcba = Mcba + mfaba * c1o3;
    mfcab = Mcab + mfaab * c1o3;
    mfbac = Mbac + mfbaa * c1o3;

    // from linearized orthogonalization 5th order central moments to central moments
    mfcbc = Mcbc + mfaba * c1o9;
    mfbcc = Mbcc + mfbaa * c1o9;
    mfccb = Mccb + mfaab * c1o9;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from  central moments to distributions as defined in Eq. (88)-(96) in \ref
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
    backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
    backwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
    backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
    backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
    backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
    backwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
    backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
    backwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
    backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
    backwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
    backwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
    backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
    backwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
    backwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
    backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
    backwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
    backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
    backwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
    backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
    backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
    backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
    backwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
    backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
    backwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

    parameters.concentration = concentration;
}
////////////////////////////////////////////////////////////////////////////////

}
//! \}
