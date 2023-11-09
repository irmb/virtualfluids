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
//! \file HighViscosityNoSlipBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "HighViscosityNoSlipBCStrategy.h"
#include "BoundaryConditions.h"
#include "DistributionArray3D.h"

HighViscosityNoSlipBCStrategy::HighViscosityNoSlipBCStrategy()
{
    BCStrategy::type         = BCStrategy::HighViscosityNoSlipBCStrategy;
    BCStrategy::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
HighViscosityNoSlipBCStrategy::~HighViscosityNoSlipBCStrategy() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> HighViscosityNoSlipBCStrategy::clone()
{
    SPtr<BCStrategy> bc(new HighViscosityNoSlipBCStrategy());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void HighViscosityNoSlipBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void HighViscosityNoSlipBCStrategy::applyBC()
{
    real f[D3Q27System::ENDF + 1];
    real feq[D3Q27System::ENDF + 1];
    distributions->getPreCollisionDistribution(f, x1, x2, x3);
    real rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);
    calcFeqFct(feq, rho, vx1, vx2, vx3);

    for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++) {
        if (bcPtr->hasNoSlipBoundaryFlag(fDir)) {
            // quadratic bounce back
            const int invDir = D3Q27System::INVDIR[fDir];
            real q        = bcPtr->getQ(invDir);
            real fReturn =
                (f[invDir] + q * f[fDir] + q * collFactor * (feq[invDir] - f[invDir] + feq[fDir] - f[fDir])) /
                (vf::basics::constant::c1o1 + q);
            distributions->setPreCollisionDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir],
                                                          x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir],
                                                          invDir);
        }
    }
}