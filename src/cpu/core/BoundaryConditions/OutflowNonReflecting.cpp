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
#include "OutflowNonReflecting.h"

#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "DistributionArray3D.h"

OutflowNonReflecting::OutflowNonReflecting()
{
    BCStrategy::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
OutflowNonReflecting::~OutflowNonReflecting() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> OutflowNonReflecting::clone()
{
    SPtr<BCStrategy> bc(new OutflowNonReflecting());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void OutflowNonReflecting::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void OutflowNonReflecting::applyBC()
{
    using namespace vf::lbm::dir;

    using namespace D3Q27System;
 //   using namespace UbMath;
    using namespace vf::basics::constant;

    real f[ENDF + 1];
    real ftemp[ENDF + 1];

    int nx1       = x1;
    int nx2       = x2;
    int nx3       = x3;
    int direction = -1;

    // flag points in direction of fluid
    if (bcPtr->hasDensityBoundaryFlag(dP00)) {
        nx1 += 1;
        direction = dP00;
    } else if (bcPtr->hasDensityBoundaryFlag(dM00)) {
        nx1 -= 1;
        direction = dM00;
    } else if (bcPtr->hasDensityBoundaryFlag(d0P0)) {
        nx2 += 1;
        direction = d0P0;
    } else if (bcPtr->hasDensityBoundaryFlag(d0M0)) {
        nx2 -= 1;
        direction = d0M0;
    } else if (bcPtr->hasDensityBoundaryFlag(d00P)) {
        nx3 += 1;
        direction = d00P;
    } else if (bcPtr->hasDensityBoundaryFlag(d00M)) {
        nx3 -= 1;
        direction = d00M;
    } else
        UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

    distributions->getPreCollisionDistribution(f, x1, x2, x3);
    distributions->getPreCollisionDistribution(ftemp, nx1, nx2, nx3);

    real rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);

    switch (direction) {
        case dP00:
            f[dP00]   = ftemp[dP00] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dP00];
            f[dPP0]  = ftemp[dPP0] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPP0];
            f[dPM0]  = ftemp[dPM0] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPM0];
            f[dP0P]  = ftemp[dP0P] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dP0P];
            f[dP0M]  = ftemp[dP0M] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dP0M];
            f[dPPP] = ftemp[dPPP] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPPP];
            f[dPMP] = ftemp[dPMP] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPMP];
            f[dPPM] = ftemp[dPPM] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPPM];
            f[dPMM] = ftemp[dPMM] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPMM];

            distributions->setPreCollisionDistributionForDirection(f[dP00], x1 + DX1[dM00], x2 + DX2[dM00], x3 + DX3[dM00], dM00);
            distributions->setPreCollisionDistributionForDirection(f[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributions->setPreCollisionDistributionForDirection(f[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributions->setPreCollisionDistributionForDirection(f[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributions->setPreCollisionDistributionForDirection(f[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            break;
        case dM00:
            f[dM00]   = ftemp[dM00] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dM00];
            f[dMP0]  = ftemp[dMP0] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMP0];
            f[dMM0]  = ftemp[dMM0] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMM0];
            f[dM0P]  = ftemp[dM0P] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dM0P];
            f[dM0M]  = ftemp[dM0M] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dM0M];
            f[dMPP] = ftemp[dMPP] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMPP];
            f[dMMP] = ftemp[dMMP] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMMP];
            f[dMPM] = ftemp[dMPM] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMPM];
            f[dMMM] = ftemp[dMMM] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[dM00], x1 + DX1[dP00], x2 + DX2[dP00], x3 + DX3[dP00], dP00);
            distributions->setPreCollisionDistributionForDirection(f[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributions->setPreCollisionDistributionForDirection(f[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributions->setPreCollisionDistributionForDirection(f[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributions->setPreCollisionDistributionForDirection(f[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);
            break;
        case d0P0:
            f[d0P0]   = ftemp[d0P0] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[d0P0];
            f[dPP0]  = ftemp[dPP0] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dPP0];
            f[dMP0]  = ftemp[dMP0] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dMP0];
            f[d0PP]  = ftemp[d0PP] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[d0PP];
            f[d0PM]  = ftemp[d0PM] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[d0PM];
            f[dPPP] = ftemp[dPPP] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dPPP];
            f[dMPP] = ftemp[dMPP] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dMPP];
            f[dPPM] = ftemp[dPPM] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dPPM];
            f[dMPM] = ftemp[dMPM] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dMPM];

            distributions->setPreCollisionDistributionForDirection(f[d0P0], x1 + DX1[d0M0], x2 + DX2[d0M0], x3 + DX3[d0M0], d0M0);
            distributions->setPreCollisionDistributionForDirection(f[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributions->setPreCollisionDistributionForDirection(f[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributions->setPreCollisionDistributionForDirection(f[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributions->setPreCollisionDistributionForDirection(f[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            break;
        case d0M0:
            f[d0M0]   = ftemp[d0M0] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[d0M0];
            f[dPM0]  = ftemp[dPM0] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dPM0];
            f[dMM0]  = ftemp[dMM0] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dMM0];
            f[d0MP]  = ftemp[d0MP] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[d0MP];
            f[d0MM]  = ftemp[d0MM] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[d0MM];
            f[dPMP] = ftemp[dPMP] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dPMP];
            f[dMMP] = ftemp[dMMP] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dMMP];
            f[dPMM] = ftemp[dPMM] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dPMM];
            f[dMMM] = ftemp[dMMM] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[d0M0], x1 + DX1[d0P0], x2 + DX2[d0P0], x3 + DX3[d0P0], d0P0);
            distributions->setPreCollisionDistributionForDirection(f[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributions->setPreCollisionDistributionForDirection(f[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributions->setPreCollisionDistributionForDirection(f[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributions->setPreCollisionDistributionForDirection(f[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);
            break;
        case d00P:
            f[d00P]   = ftemp[d00P] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[d00P];
            f[dP0P]  = ftemp[dP0P] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dP0P];
            f[dM0P]  = ftemp[dM0P] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dM0P];
            f[d0PP]  = ftemp[d0PP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[d0PP];
            f[d0MP]  = ftemp[d0MP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[d0MP];
            f[dPPP] = ftemp[dPPP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dPPP];
            f[dMPP] = ftemp[dMPP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dMPP];
            f[dPMP] = ftemp[dPMP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dPMP];
            f[dMMP] = ftemp[dMMP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dMMP];

            distributions->setPreCollisionDistributionForDirection(f[d00P], x1 + DX1[d00M], x2 + DX2[d00M], x3 + DX3[d00M], d00M);
            distributions->setPreCollisionDistributionForDirection(f[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributions->setPreCollisionDistributionForDirection(f[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributions->setPreCollisionDistributionForDirection(f[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributions->setPreCollisionDistributionForDirection(f[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            break;
        case d00M:
            f[d00M]   = ftemp[d00M] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[d00M];
            f[dP0M]  = ftemp[dP0M] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dP0M];
            f[dM0M]  = ftemp[dM0M] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dM0M];
            f[d0PM]  = ftemp[d0PM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[d0PM];
            f[d0MM]  = ftemp[d0MM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[d0MM];
            f[dPPM] = ftemp[dPPM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dPPM];
            f[dMPM] = ftemp[dMPM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dMPM];
            f[dPMM] = ftemp[dPMM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dPMM];
            f[dMMM] = ftemp[dMMM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[d00M], x1 + DX1[d00P], x2 + DX2[d00P], x3 + DX3[d00P], d00P);
            distributions->setPreCollisionDistributionForDirection(f[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributions->setPreCollisionDistributionForDirection(f[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributions->setPreCollisionDistributionForDirection(f[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributions->setPreCollisionDistributionForDirection(f[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);
            break;
        default:
            UB_THROW(
                UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
    }
}

//! \}
