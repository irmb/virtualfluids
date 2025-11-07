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
//! \addtogroup cpu_BoundaryConditions BoundaryConditions
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#include "NoSlipInterpolatedRelaxed.h"
#include "BoundaryConditions.h"
#include "DistributionArray3D.h"

NoSlipInterpolatedRelaxed::NoSlipInterpolatedRelaxed()
{
    BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> NoSlipInterpolatedRelaxed::clone()
{
    SPtr<BCStrategy> bc(new NoSlipInterpolatedRelaxed());
    std::static_pointer_cast<NoSlipInterpolatedRelaxed>(bc)->thirdMomentsFactor = thirdMomentsFactor;
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void NoSlipInterpolatedRelaxed::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void NoSlipInterpolatedRelaxed::applyBC()
{
    using namespace vf::basics::constant;
    using namespace d3q27_system;
    real f[ENDF + 1];
    real feq[ENDF + 1];
    distributions->getPostCollisionDistribution(f, x1, x2, x3);
    real rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);
    calcFeqFct(feq, rho, vx1, vx2, vx3);

    for (int fdir = FSTARTDIR; fdir <= FENDDIR; fdir++) {
        if (bcPtr->hasNoSlipBoundaryFlag(fdir)) {
            // quadratic bounce back
            const int invDir = INVDIR[fdir];
            real q = bcPtr->getQ(invDir);
            real fReturn = ((c1o1 - q) / (c1o1 + q)) * ((feq[invDir]-feq[fdir])+((f[invDir]+f[fdir])-(feq[invDir]+feq[fdir])*collFactor) / (c1o1-collFactor))*c1o2  + ((q / (c1o1 + q)) * (f[invDir] + f[fdir]))+((f[invDir]-f[fdir])-(feq[invDir]-feq[fdir]))*c1o2*thirdMomentsFactor;
            real fRetOld=((feq[fdir]-feq[invDir])+((f[invDir]+f[fdir])-(feq[invDir]+feq[fdir])*collFactor) / (c1o1-collFactor))*c1o2;
            real alpha = c1o2;//Relaxation parameter between 0 and 1
            fReturn = fReturn * alpha + (c1o1 - alpha) * fRetOld;
            distributions->setPostCollisionDistributionForDirection(fReturn, x1 + DX1[invDir], x2 + DX2[invDir], x3 + DX3[invDir], fdir);

        }
    }
}

void NoSlipInterpolatedRelaxed::thirdMomentsOn()
{
    this->thirdMomentsFactor = 1.0;
}

void NoSlipInterpolatedRelaxed::thirdMomentsOff()
{
    this->thirdMomentsFactor = 0.0;
}

//! \}
