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
#include "VelocityWithPressureInterpolated.h"
#include "BCArray3D.h"
#include "DistributionArray3D.h"

VelocityWithPressureInterpolated::VelocityWithPressureInterpolated()
{
    BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
VelocityWithPressureInterpolated::~VelocityWithPressureInterpolated() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> VelocityWithPressureInterpolated::clone()
{
    SPtr<BCStrategy> bc(new VelocityWithPressureInterpolated());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void VelocityWithPressureInterpolated::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void VelocityWithPressureInterpolated::applyBC()
{
    using namespace vf::basics::constant;

   //velocity bc for non reflecting pressure bc
   real f[D3Q27System::ENDF+1];
   //real feq[D3Q27System::ENDF+1];
   distributions->getPostCollisionDistribution(f, x1, x2, x3);
   real rho, vx1, vx2, vx3, drho;
   calcMacrosFct(f, drho, vx1, vx2, vx3);
   //calcFeqFct(feq, drho, vx1, vx2, vx3);
   
   rho = c1o1+drho*compressibleFactor;

   for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
   {
        int nX1 = x1 + D3Q27System::DX1[fdir];
        int nX2 = x2 + D3Q27System::DX2[fdir];
        int nX3 = x3 + D3Q27System::DX3[fdir];

        int minX1 = 0;
        int minX2 = 0;
        int minX3 = 0;

        int maxX1 = (int)bcArray->getNX1();
        int maxX2 = (int)bcArray->getNX2();
        int maxX3 = (int)bcArray->getNX3();

        if (minX1 <= nX1 && maxX1 > nX1 && minX2 <= nX2 && maxX2 > nX2 && minX3 <= nX3 && maxX3 > nX3) {
            if (bcArray->isSolid(nX1, nX2, nX3)) {
                const int invDir = D3Q27System::INVDIR[fdir];
                //real q =1.0;// bcPtr->getQ(invDir);// m+m q=0 stabiler
                real velocity = bcPtr->getBoundaryVelocity(fdir);
                
                //real fReturn = ((1.0 - q) / (1.0 + q))*((f[fdir] - feq[fdir]*collFactor) / (1.0 -
                //collFactor)) + ((q*(f[fdir] + f[invDir]) - velocity*rho) / (1.0 +
                //q))-drho*D3Q27System::WEIGTH[invDir];

                // if q=1
                // real fReturn = ((q*(f[fdir] + f[invDir]) - velocity*rho) / (1.0 +
                // q))-drho*D3Q27System::WEIGTH[invDir];
                real fReturn = (f[fdir] + f[invDir] - velocity * rho) / vf::basics::constant::c2o1 - drho * D3Q27System::WEIGTH[invDir];

                distributions->setPostCollisionDistributionForDirection(fReturn, nX1, nX2, nX3, invDir);
            }
        }
    }
}

//! \}
