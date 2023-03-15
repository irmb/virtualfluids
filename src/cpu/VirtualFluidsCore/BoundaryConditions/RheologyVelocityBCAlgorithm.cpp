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
//! \file RheologyVelocityBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "RheologyVelocityBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

RheologyVelocityBCAlgorithm::RheologyVelocityBCAlgorithm()
{
   //BCAlgorithm::type = BCAlgorithm::RheologyVelocityBCAlgorithm;
   //BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
RheologyVelocityBCAlgorithm::~RheologyVelocityBCAlgorithm()
{
}
//////////////////////////////////////////////////////////////////////////
void RheologyVelocityBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void RheologyVelocityBCAlgorithm::applyBC()
{
    using namespace vf::lbm::constant;

    real f[D3Q27System::ENDF+1];
    real feq[D3Q27System::ENDF+1];
    distributions->getDistribution(f, x1, x2, x3);
    real rho, vx1, vx2, vx3, drho;
    calcMacrosFct(f, drho, vx1, vx2, vx3);
    calcFeqFct(feq, drho, vx1, vx2, vx3);

    real shearRate = D3Q27System::getShearRate(f, collFactor);
    // getShearRate(f, collFactor) takes pre collision state but is given post collision state. Requires rescaling
    real collFactorF = getRheologyCollFactor(collFactor, shearRate, drho);


    rho = 1.0+drho*compressibleFactor;

    for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
    {
        if (bcPtr->hasVelocityBoundaryFlag(fdir))
        {
            const int invDir = D3Q27System::INVDIR[fdir];
            real q = bcPtr->getQ(invDir);
            real velocity = bcPtr->getBoundaryVelocity(invDir);
            real fReturn = ((c1o1 - q) / (c1o1 + q)) * ((f[invDir])) + ((q * (f[invDir] + f[fdir] - collFactorF * (f[invDir] - feq[invDir] + f[fdir] - feq[fdir])) - velocity * rho) / (c1o1 + q));
            distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
        }
    }

}

