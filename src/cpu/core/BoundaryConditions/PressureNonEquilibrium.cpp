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
#include "PressureNonEquilibrium.h"

#include "BoundaryConditions.h"
#include "DistributionArray3D.h"

PressureNonEquilibrium::PressureNonEquilibrium()
{
    BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
PressureNonEquilibrium::~PressureNonEquilibrium() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> PressureNonEquilibrium::clone()
{
    SPtr<BCStrategy> bc(new PressureNonEquilibrium());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void PressureNonEquilibrium::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void PressureNonEquilibrium::applyBC()
{
    using namespace vf::lbm::dir;

    real f[D3Q27System::ENDF + 1];
    distributions->getPostCollisionDistribution(f, x1, x2, x3);
    int nx1 = x1;
    int nx2 = x2;
    int nx3 = x3;

    // flag points in direction of fluid
    if (bcPtr->hasDensityBoundaryFlag(dP00)) {
        nx1 -= 1;
    } else if (bcPtr->hasDensityBoundaryFlag(dM00)) {
        nx1 += 1;
    } else if (bcPtr->hasDensityBoundaryFlag(d0P0)) {
        nx2 -= 1;
    } else if (bcPtr->hasDensityBoundaryFlag(d0M0)) {
        nx2 += 1;
    } else if (bcPtr->hasDensityBoundaryFlag(d00P)) {
        nx3 -= 1;
    } else if (bcPtr->hasDensityBoundaryFlag(d00M)) {
        nx3 += 1;
    } else
        return; // UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

    real rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);
    // real vlimit=0.01;
    // vx1=(fabs(vx1)>vlimit) ? vx1/fabs(vx1)*vlimit : vx1;
    // vx2=(fabs(vx2)>vlimit) ? vx2/fabs(vx2)*vlimit : vx2;
    // vx3=(fabs(vx3)>vlimit) ? vx3/fabs(vx3)*vlimit : vx3;
    real rhoBC = bcPtr->getBoundaryDensity();
    for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++) {
        if (bcPtr->hasDensityBoundaryFlag(fdir)) {
            // Martins NEQ ADDON
            ////original: 15.2.2013:
            real ftemp = calcFeqsForDirFct(fdir, rho, vx1, vx2, vx3);
            // rhoBC=(rho>rhoBC)? rhoBC : rho; //Limiter 08.08.2018
            ftemp = calcFeqsForDirFct(fdir, rhoBC, vx1, vx2, vx3) + f[fdir] - ftemp;
            distributions->setPostCollisionDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);
        }
    }
}

//! \}
