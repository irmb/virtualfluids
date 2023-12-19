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

#include "VelocityInterpolated.h"
#include "BoundaryConditions.h"
#include "DistributionArray3D.h"
#include "Block3D.h"

VelocityInterpolated::VelocityInterpolated()
{
    BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> VelocityInterpolated::clone()
{
    SPtr<BCStrategy> bc(new VelocityInterpolated());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void VelocityInterpolated::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void VelocityInterpolated::applyBC()
{
    real f[D3Q27System::ENDF + 1];
    real feq[D3Q27System::ENDF + 1];
    distributions->getPostCollisionDistribution(f, x1, x2, x3);
    real rho, vx1, vx2, vx3, drho;
    calcMacrosFct(f, drho, vx1, vx2, vx3);
    calcFeqFct(feq, drho, vx1, vx2, vx3);

    //DEBUG
    //int blockID = block->getGlobalID();

    rho = vf::basics::constant::c1o1 + drho * compressibleFactor;

    for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
        if (bcPtr->hasVelocityBoundaryFlag(fdir)) {
            const int invDir = D3Q27System::INVDIR[fdir];
            real q        = bcPtr->getQ(invDir);
            real velocity = bcPtr->getBoundaryVelocity(invDir);
            real fReturn = ((vf::basics::constant::c1o1 - q) / (vf::basics::constant::c1o1 + q)) * ((f[invDir] - feq[invDir]) / (vf::basics::constant::c1o1 - collFactor) + feq[invDir]) +
                              ((q * (f[invDir] + f[fdir]) - velocity * rho) / (vf::basics::constant::c1o1 + q));
            distributions->setPostCollisionDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir],
                                                       x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir],
                                                       fdir);
        }
    }
}

//! \}
