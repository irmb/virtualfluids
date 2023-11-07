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
//! \file NonReflectingOutflowBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "NonReflectingOutflowBCStrategy.h"

#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "DistributionArray3D.h"

NonReflectingOutflowBCStrategy::NonReflectingOutflowBCStrategy()
{
    BCStrategy::type         = BCStrategy::NonReflectingOutflowBCStrategy;
    BCStrategy::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
NonReflectingOutflowBCStrategy::~NonReflectingOutflowBCStrategy() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> NonReflectingOutflowBCStrategy::clone()
{
    SPtr<BCStrategy> bc(new NonReflectingOutflowBCStrategy());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void NonReflectingOutflowBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void NonReflectingOutflowBCStrategy::applyBC()
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
    } else if (bcPtr->hasDensityBoundaryFlag(DIR_0P0)) {
        nx2 += 1;
        direction = DIR_0P0;
    } else if (bcPtr->hasDensityBoundaryFlag(DIR_0M0)) {
        nx2 -= 1;
        direction = DIR_0M0;
    } else if (bcPtr->hasDensityBoundaryFlag(DIR_00P)) {
        nx3 += 1;
        direction = DIR_00P;
    } else if (bcPtr->hasDensityBoundaryFlag(DIR_00M)) {
        nx3 -= 1;
        direction = DIR_00M;
    } else
        UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

    distributions->getPreCollisionDistribution(f, x1, x2, x3);
    distributions->getPreCollisionDistribution(ftemp, nx1, nx2, nx3);

    real rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);

    switch (direction) {
        case dP00:
            f[dP00]   = ftemp[dP00] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dP00];
            f[DIR_PP0]  = ftemp[DIR_PP0] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[DIR_PP0];
            f[DIR_PM0]  = ftemp[DIR_PM0] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[DIR_PM0];
            f[DIR_P0P]  = ftemp[DIR_P0P] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[DIR_P0P];
            f[DIR_P0M]  = ftemp[DIR_P0M] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[DIR_P0M];
            f[DIR_PPP] = ftemp[DIR_PPP] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[DIR_PPP];
            f[DIR_PMP] = ftemp[DIR_PMP] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[DIR_PMP];
            f[DIR_PPM] = ftemp[DIR_PPM] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[DIR_PPM];
            f[DIR_PMM] = ftemp[DIR_PMM] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[DIR_PMM];

            distributions->setPreCollisionDistributionForDirection(f[dP00], x1 + DX1[dM00], x2 + DX2[dM00], x3 + DX3[dM00], dM00);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PP0], x1 + DX1[DIR_MM0], x2 + DX2[DIR_MM0], x3 + DX3[DIR_MM0], DIR_MM0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PM0], x1 + DX1[DIR_MP0], x2 + DX2[DIR_MP0], x3 + DX3[DIR_MP0], DIR_MP0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_P0P], x1 + DX1[DIR_M0M], x2 + DX2[DIR_M0M], x3 + DX3[DIR_M0M], DIR_M0M);
            distributions->setPreCollisionDistributionForDirection(f[DIR_P0M], x1 + DX1[DIR_M0P], x2 + DX2[DIR_M0P], x3 + DX3[DIR_M0P], DIR_M0P);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
            break;
        case dM00:
            f[dM00]   = ftemp[dM00] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dM00];
            f[DIR_MP0]  = ftemp[DIR_MP0] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[DIR_MP0];
            f[DIR_MM0]  = ftemp[DIR_MM0] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[DIR_MM0];
            f[DIR_M0P]  = ftemp[DIR_M0P] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[DIR_M0P];
            f[DIR_M0M]  = ftemp[DIR_M0M] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[DIR_M0M];
            f[DIR_MPP] = ftemp[DIR_MPP] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[DIR_MPP];
            f[DIR_MMP] = ftemp[DIR_MMP] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[DIR_MMP];
            f[DIR_MPM] = ftemp[DIR_MPM] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[DIR_MPM];
            f[DIR_MMM] = ftemp[DIR_MMM] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[DIR_MMM];

            distributions->setPreCollisionDistributionForDirection(f[dM00], x1 + DX1[dP00], x2 + DX2[dP00], x3 + DX3[dP00], dP00);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MP0], x1 + DX1[DIR_PM0], x2 + DX2[DIR_PM0], x3 + DX3[DIR_PM0], DIR_PM0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MM0], x1 + DX1[DIR_PP0], x2 + DX2[DIR_PP0], x3 + DX3[DIR_PP0], DIR_PP0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_M0P], x1 + DX1[DIR_P0M], x2 + DX2[DIR_P0M], x3 + DX3[DIR_P0M], DIR_P0M);
            distributions->setPreCollisionDistributionForDirection(f[DIR_M0M], x1 + DX1[DIR_P0P], x2 + DX2[DIR_P0P], x3 + DX3[DIR_P0P], DIR_P0P);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
            break;
        case DIR_0P0:
            f[DIR_0P0]   = ftemp[DIR_0P0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_0P0];
            f[DIR_PP0]  = ftemp[DIR_PP0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_PP0];
            f[DIR_MP0]  = ftemp[DIR_MP0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_MP0];
            f[DIR_0PP]  = ftemp[DIR_0PP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_0PP];
            f[DIR_0PM]  = ftemp[DIR_0PM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_0PM];
            f[DIR_PPP] = ftemp[DIR_PPP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_PPP];
            f[DIR_MPP] = ftemp[DIR_MPP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_MPP];
            f[DIR_PPM] = ftemp[DIR_PPM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_PPM];
            f[DIR_MPM] = ftemp[DIR_MPM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[DIR_MPM];

            distributions->setPreCollisionDistributionForDirection(f[DIR_0P0], x1 + DX1[DIR_0M0], x2 + DX2[DIR_0M0], x3 + DX3[DIR_0M0], DIR_0M0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PP0], x1 + DX1[DIR_MM0], x2 + DX2[DIR_MM0], x3 + DX3[DIR_MM0], DIR_MM0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MP0], x1 + DX1[DIR_PM0], x2 + DX2[DIR_PM0], x3 + DX3[DIR_PM0], DIR_PM0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_0PP], x1 + DX1[DIR_0MM], x2 + DX2[DIR_0MM], x3 + DX3[DIR_0MM], DIR_0MM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_0PM], x1 + DX1[DIR_0MP], x2 + DX2[DIR_0MP], x3 + DX3[DIR_0MP], DIR_0MP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
            break;
        case DIR_0M0:
            f[DIR_0M0]   = ftemp[DIR_0M0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_0M0];
            f[DIR_PM0]  = ftemp[DIR_PM0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_PM0];
            f[DIR_MM0]  = ftemp[DIR_MM0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_MM0];
            f[DIR_0MP]  = ftemp[DIR_0MP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_0MP];
            f[DIR_0MM]  = ftemp[DIR_0MM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_0MM];
            f[DIR_PMP] = ftemp[DIR_PMP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_PMP];
            f[DIR_MMP] = ftemp[DIR_MMP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_MMP];
            f[DIR_PMM] = ftemp[DIR_PMM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_PMM];
            f[DIR_MMM] = ftemp[DIR_MMM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[DIR_MMM];

            distributions->setPreCollisionDistributionForDirection(f[DIR_0M0], x1 + DX1[DIR_0P0], x2 + DX2[DIR_0P0], x3 + DX3[DIR_0P0], DIR_0P0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PM0], x1 + DX1[DIR_MP0], x2 + DX2[DIR_MP0], x3 + DX3[DIR_MP0], DIR_MP0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MM0], x1 + DX1[DIR_PP0], x2 + DX2[DIR_PP0], x3 + DX3[DIR_PP0], DIR_PP0);
            distributions->setPreCollisionDistributionForDirection(f[DIR_0MP], x1 + DX1[DIR_0PM], x2 + DX2[DIR_0PM], x3 + DX3[DIR_0PM], DIR_0PM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_0MM], x1 + DX1[DIR_0PP], x2 + DX2[DIR_0PP], x3 + DX3[DIR_0PP], DIR_0PP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
            break;
        case DIR_00P:
            f[DIR_00P]   = ftemp[DIR_00P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_00P];
            f[DIR_P0P]  = ftemp[DIR_P0P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_P0P];
            f[DIR_M0P]  = ftemp[DIR_M0P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_M0P];
            f[DIR_0PP]  = ftemp[DIR_0PP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_0PP];
            f[DIR_0MP]  = ftemp[DIR_0MP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_0MP];
            f[DIR_PPP] = ftemp[DIR_PPP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_PPP];
            f[DIR_MPP] = ftemp[DIR_MPP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_MPP];
            f[DIR_PMP] = ftemp[DIR_PMP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_PMP];
            f[DIR_MMP] = ftemp[DIR_MMP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[DIR_MMP];

            distributions->setPreCollisionDistributionForDirection(f[DIR_00P], x1 + DX1[DIR_00M], x2 + DX2[DIR_00M], x3 + DX3[DIR_00M], DIR_00M);
            distributions->setPreCollisionDistributionForDirection(f[DIR_P0P], x1 + DX1[DIR_M0M], x2 + DX2[DIR_M0M], x3 + DX3[DIR_M0M], DIR_M0M);
            distributions->setPreCollisionDistributionForDirection(f[DIR_M0P], x1 + DX1[DIR_P0M], x2 + DX2[DIR_P0M], x3 + DX3[DIR_P0M], DIR_P0M);
            distributions->setPreCollisionDistributionForDirection(f[DIR_0PP], x1 + DX1[DIR_0MM], x2 + DX2[DIR_0MM], x3 + DX3[DIR_0MM], DIR_0MM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_0MP], x1 + DX1[DIR_0PM], x2 + DX2[DIR_0PM], x3 + DX3[DIR_0PM], DIR_0PM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
            break;
        case DIR_00M:
            f[DIR_00M]   = ftemp[DIR_00M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_00M];
            f[DIR_P0M]  = ftemp[DIR_P0M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_P0M];
            f[DIR_M0M]  = ftemp[DIR_M0M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_M0M];
            f[DIR_0PM]  = ftemp[DIR_0PM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_0PM];
            f[DIR_0MM]  = ftemp[DIR_0MM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_0MM];
            f[DIR_PPM] = ftemp[DIR_PPM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_PPM];
            f[DIR_MPM] = ftemp[DIR_MPM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_MPM];
            f[DIR_PMM] = ftemp[DIR_PMM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_PMM];
            f[DIR_MMM] = ftemp[DIR_MMM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[DIR_MMM];

            distributions->setPreCollisionDistributionForDirection(f[DIR_00M], x1 + DX1[DIR_00P], x2 + DX2[DIR_00P], x3 + DX3[DIR_00P], DIR_00P);
            distributions->setPreCollisionDistributionForDirection(f[DIR_P0M], x1 + DX1[DIR_M0P], x2 + DX2[DIR_M0P], x3 + DX3[DIR_M0P], DIR_M0P);
            distributions->setPreCollisionDistributionForDirection(f[DIR_M0M], x1 + DX1[DIR_P0P], x2 + DX2[DIR_P0P], x3 + DX3[DIR_P0P], DIR_P0P);
            distributions->setPreCollisionDistributionForDirection(f[DIR_0PM], x1 + DX1[DIR_0MP], x2 + DX2[DIR_0MP], x3 + DX3[DIR_0MP], DIR_0MP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_0MM], x1 + DX1[DIR_0PP], x2 + DX2[DIR_0PP], x3 + DX3[DIR_0PP], DIR_0PP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
            distributions->setPreCollisionDistributionForDirection(f[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
            break;
        default:
            UB_THROW(
                UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
    }
}
