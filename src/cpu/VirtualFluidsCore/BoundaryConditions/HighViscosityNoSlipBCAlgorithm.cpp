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
//! \file HighViscosityNoSlipBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "HighViscosityNoSlipBCAlgorithm.h"
#include "BoundaryConditions.h"
#include "DistributionArray3D.h"

HighViscosityNoSlipBCAlgorithm::HighViscosityNoSlipBCAlgorithm()
{
    BCAlgorithm::type         = BCAlgorithm::HighViscosityNoSlipBCAlgorithm;
    BCAlgorithm::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
HighViscosityNoSlipBCAlgorithm::~HighViscosityNoSlipBCAlgorithm() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> HighViscosityNoSlipBCAlgorithm::clone()
{
    SPtr<BCAlgorithm> bc(new HighViscosityNoSlipBCAlgorithm());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void HighViscosityNoSlipBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void HighViscosityNoSlipBCAlgorithm::applyBC()
{
    LBMReal f[D3Q27System::ENDF + 1];
    LBMReal feq[D3Q27System::ENDF + 1];
    distributions->getDistribution(f, x1, x2, x3);
    LBMReal rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);
    calcFeqFct(feq, rho, vx1, vx2, vx3);

    for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++) {
        if (bcPtr->hasNoSlipBoundaryFlag(fDir)) {
            // quadratic bounce back
            const int invDir = D3Q27System::INVDIR[fDir];
            LBMReal q        = bcPtr->getQ(invDir);
            LBMReal fReturn =
                (f[invDir] + q * f[fDir] + q * collFactor * (feq[invDir] - f[invDir] + feq[fDir] - f[fDir])) /
                (1.0 + q);
            distributions->setDistributionInvForDirection(fReturn, x1 + D3Q27System::DX1[invDir],
                                                          x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir],
                                                          invDir);
        }
    }
}
