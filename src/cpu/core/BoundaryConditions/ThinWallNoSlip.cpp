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
#include "ThinWallNoSlip.h"

#include "BoundaryConditions.h"
#include "EsoSplit.h"

ThinWallNoSlip::ThinWallNoSlip()
{
    BCStrategy::preCollision = false;
    pass                      = 1;
}
//////////////////////////////////////////////////////////////////////////
ThinWallNoSlip::~ThinWallNoSlip() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> ThinWallNoSlip::clone()
{
    SPtr<BCStrategy> bc(new ThinWallNoSlip());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlip::applyBC()
{
    real f[D3Q27System::ENDF + 1];
    real feq[D3Q27System::ENDF + 1];
    distributions->getPostCollisionDistribution(f, x1, x2, x3);
    real rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);
    calcFeqFct(feq, rho, vx1, vx2, vx3);

    real fReturn;

    for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
        if (bcPtr->hasNoSlipBoundaryFlag(fdir)) {
            const int invDir = D3Q27System::INVDIR[fdir];
            if (pass == 1) {
                real q = bcPtr->getQ(invDir);
                fReturn   = ((vf::basics::constant::c1o1 - q) / (vf::basics::constant::c1o1 + q)) * vf::basics::constant::c1o2 *
                          (f[invDir] - f[fdir] +
                           (f[invDir] + f[fdir] - collFactor * (feq[fdir] + feq[invDir])) / (vf::basics::constant::c1o1 - collFactor));
                // distributionsTemp->setPostCollisionDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 +
                // D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
                fTemp[fdir] = fReturn;
            } else {
                // quadratic bounce back with for thin walls
                // fReturn = distributionsTemp->getPostCollisionDistributionForDirection(x1 + D3Q27System::DX1[invDir], x2 +
                // D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
                fReturn = fTemp[fdir];
                distributions->setPostCollisionDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir],
                                                           x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir],
                                                           fdir);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlip::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThinWallNoSlip::setPass(int pass)
{
    this->pass = pass;
}

bool ThinWallNoSlip::isThinWallNoSlipBCStrategy()
{
    return true;
}

//! \}
