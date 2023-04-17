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
//! \file NonReflectingOutflowBCAlgorithmWithRelaxation.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher, Hussein Alihussein
//=======================================================================================
#include "NonReflectingOutflowBCAlgorithmWithRelaxation.h"

#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "DistributionArray3D.h"

NonReflectingOutflowBCAlgorithmWithRelaxation::NonReflectingOutflowBCAlgorithmWithRelaxation()
{
    BCAlgorithm::type         = BCAlgorithm::NonReflectingOutflowBCAlgorithmWithRelaxation;
    BCAlgorithm::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
NonReflectingOutflowBCAlgorithmWithRelaxation::~NonReflectingOutflowBCAlgorithmWithRelaxation() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> NonReflectingOutflowBCAlgorithmWithRelaxation::clone()
{
    SPtr<BCAlgorithm> bc(new NonReflectingOutflowBCAlgorithmWithRelaxation());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void NonReflectingOutflowBCAlgorithmWithRelaxation::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void NonReflectingOutflowBCAlgorithmWithRelaxation::applyBC()
{
    using namespace vf::lbm::dir;

    using namespace D3Q27System;
 //   using namespace UbMath;
    using namespace vf::basics::constant;

    LBMReal f[ENDF + 1];
    LBMReal ftemp[ENDF + 1];

    int nx1       = x1;
    int nx2       = x2;
    int nx3       = x3;
    int direction = -1;

    // flag points in direction of fluid
    if (bcPtr->hasDensityBoundaryFlag(DIR_P00)) {
        nx1 += 1;
        direction = DIR_P00;
    } else if (bcPtr->hasDensityBoundaryFlag(DIR_M00)) {
        nx1 -= 1;
        direction = DIR_M00;
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

    distributions->getDistribution(f, x1, x2, x3);
    distributions->getDistribution(ftemp, nx1, nx2, nx3);

    LBMReal rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);
    LBMReal delf = rho*0.01;
    switch (direction) {
        case DIR_P00:
            f[DIR_P00]   = ftemp[DIR_P00] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_P00] - delf* WEIGTH[DIR_P00];
            f[DIR_PP0]  = ftemp[DIR_PP0] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_PP0]- delf* WEIGTH[DIR_PP0];
            f[DIR_PM0]  = ftemp[DIR_PM0] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_PM0]- delf* WEIGTH[DIR_PM0];
            f[DIR_P0P]  = ftemp[DIR_P0P] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_P0P]- delf* WEIGTH[DIR_P0P];
            f[DIR_P0M]  = ftemp[DIR_P0M] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_P0M]- delf* WEIGTH[DIR_P0M];
            f[DIR_PPP] = ftemp[DIR_PPP] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_PPP]- delf* WEIGTH[DIR_PPP];
            f[DIR_PMP] = ftemp[DIR_PMP] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_PMP]- delf* WEIGTH[DIR_PMP];
            f[DIR_PPM] = ftemp[DIR_PPM] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_PPM]- delf* WEIGTH[DIR_PPM];
            f[DIR_PMM] = ftemp[DIR_PMM] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1) * f[DIR_PMM]- delf* WEIGTH[DIR_PMM];

            distributions->setDistributionInvForDirection(f[DIR_P00], x1 + DX1[DIR_M00], x2 + DX2[DIR_M00], x3 + DX3[DIR_M00], DIR_M00);
            distributions->setDistributionInvForDirection(f[DIR_PP0], x1 + DX1[DIR_MM0], x2 + DX2[DIR_MM0], x3 + DX3[DIR_MM0], DIR_MM0);
            distributions->setDistributionInvForDirection(f[DIR_PM0], x1 + DX1[DIR_MP0], x2 + DX2[DIR_MP0], x3 + DX3[DIR_MP0], DIR_MP0);
            distributions->setDistributionInvForDirection(f[DIR_P0P], x1 + DX1[DIR_M0M], x2 + DX2[DIR_M0M], x3 + DX3[DIR_M0M], DIR_M0M);
            distributions->setDistributionInvForDirection(f[DIR_P0M], x1 + DX1[DIR_M0P], x2 + DX2[DIR_M0P], x3 + DX3[DIR_M0P], DIR_M0P);
            distributions->setDistributionInvForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
            distributions->setDistributionInvForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
            distributions->setDistributionInvForDirection(f[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
            distributions->setDistributionInvForDirection(f[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
            break;
        case DIR_M00:
            f[DIR_M00]   = ftemp[DIR_M00] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_M00]- delf* WEIGTH[DIR_M00];
            f[DIR_MP0]  = ftemp[DIR_MP0] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_MP0]- delf* WEIGTH[DIR_MP0];
            f[DIR_MM0]  = ftemp[DIR_MM0] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_MM0]- delf* WEIGTH[DIR_MM0];
            f[DIR_M0P]  = ftemp[DIR_M0P] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_M0P]- delf* WEIGTH[DIR_M0P];
            f[DIR_M0M]  = ftemp[DIR_M0M] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_M0M]- delf* WEIGTH[DIR_M0M];
            f[DIR_MPP] = ftemp[DIR_MPP] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_MPP]- delf* WEIGTH[DIR_MPP];
            f[DIR_MMP] = ftemp[DIR_MMP] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_MMP]- delf* WEIGTH[DIR_MMP];
            f[DIR_MPM] = ftemp[DIR_MPM] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_MPM]- delf* WEIGTH[DIR_MPM];
            f[DIR_MMM] = ftemp[DIR_MMM] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1) * f[DIR_MMM]- delf* WEIGTH[DIR_MMM];

            distributions->setDistributionInvForDirection(f[DIR_M00], x1 + DX1[DIR_P00], x2 + DX2[DIR_P00], x3 + DX3[DIR_P00], DIR_P00);
            distributions->setDistributionInvForDirection(f[DIR_MP0], x1 + DX1[DIR_PM0], x2 + DX2[DIR_PM0], x3 + DX3[DIR_PM0], DIR_PM0);
            distributions->setDistributionInvForDirection(f[DIR_MM0], x1 + DX1[DIR_PP0], x2 + DX2[DIR_PP0], x3 + DX3[DIR_PP0], DIR_PP0);
            distributions->setDistributionInvForDirection(f[DIR_M0P], x1 + DX1[DIR_P0M], x2 + DX2[DIR_P0M], x3 + DX3[DIR_P0M], DIR_P0M);
            distributions->setDistributionInvForDirection(f[DIR_M0M], x1 + DX1[DIR_P0P], x2 + DX2[DIR_P0P], x3 + DX3[DIR_P0P], DIR_P0P);
            distributions->setDistributionInvForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
            distributions->setDistributionInvForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
            distributions->setDistributionInvForDirection(f[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
            distributions->setDistributionInvForDirection(f[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
            break;
        case DIR_0P0:
            f[DIR_0P0]   = ftemp[DIR_0P0] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_0P0]- delf* WEIGTH[DIR_0P0];
            f[DIR_PP0]  = ftemp[DIR_PP0] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_PP0]- delf* WEIGTH[DIR_PP0];
            f[DIR_MP0]  = ftemp[DIR_MP0] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_MP0]- delf* WEIGTH[DIR_MP0];
            f[DIR_0PP]  = ftemp[DIR_0PP] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_0PP]- delf* WEIGTH[DIR_0PP];
            f[DIR_0PM]  = ftemp[DIR_0PM] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_0PM]- delf* WEIGTH[DIR_0PM];
            f[DIR_PPP] = ftemp[DIR_PPP] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_PPP]- delf* WEIGTH[DIR_PPP];
            f[DIR_MPP] = ftemp[DIR_MPP] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_MPP]- delf* WEIGTH[DIR_MPP];
            f[DIR_PPM] = ftemp[DIR_PPM] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_PPM]- delf* WEIGTH[DIR_PPM];
            f[DIR_MPM] = ftemp[DIR_MPM] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2) * f[DIR_MPM]- delf* WEIGTH[DIR_MPM];

            distributions->setDistributionInvForDirection(f[DIR_0P0], x1 + DX1[DIR_0M0], x2 + DX2[DIR_0M0], x3 + DX3[DIR_0M0], DIR_0M0);
            distributions->setDistributionInvForDirection(f[DIR_PP0], x1 + DX1[DIR_MM0], x2 + DX2[DIR_MM0], x3 + DX3[DIR_MM0], DIR_MM0);
            distributions->setDistributionInvForDirection(f[DIR_MP0], x1 + DX1[DIR_PM0], x2 + DX2[DIR_PM0], x3 + DX3[DIR_PM0], DIR_PM0);
            distributions->setDistributionInvForDirection(f[DIR_0PP], x1 + DX1[DIR_0MM], x2 + DX2[DIR_0MM], x3 + DX3[DIR_0MM], DIR_0MM);
            distributions->setDistributionInvForDirection(f[DIR_0PM], x1 + DX1[DIR_0MP], x2 + DX2[DIR_0MP], x3 + DX3[DIR_0MP], DIR_0MP);
            distributions->setDistributionInvForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
            distributions->setDistributionInvForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
            distributions->setDistributionInvForDirection(f[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
            distributions->setDistributionInvForDirection(f[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
            break;
        case DIR_0M0:
            f[DIR_0M0]   = ftemp[DIR_0M0] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_0M0]- delf* WEIGTH[DIR_0M0];
            f[DIR_PM0]  = ftemp[DIR_PM0] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_PM0]- delf* WEIGTH[DIR_PM0];
            f[DIR_MM0]  = ftemp[DIR_MM0] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_MM0]- delf* WEIGTH[DIR_MM0];
            f[DIR_0MP]  = ftemp[DIR_0MP] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_0MP]- delf* WEIGTH[DIR_0MP];
            f[DIR_0MM]  = ftemp[DIR_0MM] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_0MM]- delf* WEIGTH[DIR_0MM];
            f[DIR_PMP] = ftemp[DIR_PMP] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_PMP]- delf* WEIGTH[DIR_PMP];
            f[DIR_MMP] = ftemp[DIR_MMP] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_MMP]- delf* WEIGTH[DIR_MMP];
            f[DIR_PMM] = ftemp[DIR_PMM] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_PMM]- delf* WEIGTH[DIR_PMM];
            f[DIR_MMM] = ftemp[DIR_MMM] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2) * f[DIR_MMM]- delf* WEIGTH[DIR_MMM];

            distributions->setDistributionInvForDirection(f[DIR_0M0], x1 + DX1[DIR_0P0], x2 + DX2[DIR_0P0], x3 + DX3[DIR_0P0], DIR_0P0);
            distributions->setDistributionInvForDirection(f[DIR_PM0], x1 + DX1[DIR_MP0], x2 + DX2[DIR_MP0], x3 + DX3[DIR_MP0], DIR_MP0);
            distributions->setDistributionInvForDirection(f[DIR_MM0], x1 + DX1[DIR_PP0], x2 + DX2[DIR_PP0], x3 + DX3[DIR_PP0], DIR_PP0);
            distributions->setDistributionInvForDirection(f[DIR_0MP], x1 + DX1[DIR_0PM], x2 + DX2[DIR_0PM], x3 + DX3[DIR_0PM], DIR_0PM);
            distributions->setDistributionInvForDirection(f[DIR_0MM], x1 + DX1[DIR_0PP], x2 + DX2[DIR_0PP], x3 + DX3[DIR_0PP], DIR_0PP);
            distributions->setDistributionInvForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
            distributions->setDistributionInvForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
            distributions->setDistributionInvForDirection(f[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
            distributions->setDistributionInvForDirection(f[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
            break;
        case DIR_00P:
            f[DIR_00P]   = ftemp[DIR_00P] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_00P]- delf* WEIGTH[DIR_00P];
            f[DIR_P0P]  = ftemp[DIR_P0P] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_P0P]- delf* WEIGTH[DIR_P0P];
            f[DIR_M0P]  = ftemp[DIR_M0P] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_M0P]- delf* WEIGTH[DIR_M0P];
            f[DIR_0PP]  = ftemp[DIR_0PP] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_0PP]- delf* WEIGTH[DIR_0PP];
            f[DIR_0MP]  = ftemp[DIR_0MP] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_0MP]- delf* WEIGTH[DIR_0MP];
            f[DIR_PPP] = ftemp[DIR_PPP] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_PPP]- delf* WEIGTH[DIR_PPP];
            f[DIR_MPP] = ftemp[DIR_MPP] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_MPP]- delf* WEIGTH[DIR_MPP];
            f[DIR_PMP] = ftemp[DIR_PMP] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_PMP]- delf* WEIGTH[DIR_PMP];
            f[DIR_MMP] = ftemp[DIR_MMP] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3) * f[DIR_MMP]- delf* WEIGTH[DIR_MMP];

            distributions->setDistributionInvForDirection(f[DIR_00P], x1 + DX1[DIR_00M], x2 + DX2[DIR_00M], x3 + DX3[DIR_00M], DIR_00M);
            distributions->setDistributionInvForDirection(f[DIR_P0P], x1 + DX1[DIR_M0M], x2 + DX2[DIR_M0M], x3 + DX3[DIR_M0M], DIR_M0M);
            distributions->setDistributionInvForDirection(f[DIR_M0P], x1 + DX1[DIR_P0M], x2 + DX2[DIR_P0M], x3 + DX3[DIR_P0M], DIR_P0M);
            distributions->setDistributionInvForDirection(f[DIR_0PP], x1 + DX1[DIR_0MM], x2 + DX2[DIR_0MM], x3 + DX3[DIR_0MM], DIR_0MM);
            distributions->setDistributionInvForDirection(f[DIR_0MP], x1 + DX1[DIR_0PM], x2 + DX2[DIR_0PM], x3 + DX3[DIR_0PM], DIR_0PM);
            distributions->setDistributionInvForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
            distributions->setDistributionInvForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
            distributions->setDistributionInvForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
            distributions->setDistributionInvForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
            break;
        case DIR_00M:
            f[DIR_00M]   = ftemp[DIR_00M] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_00M]- delf* WEIGTH[DIR_00M];
            f[DIR_P0M]  = ftemp[DIR_P0M] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_P0M]- delf* WEIGTH[DIR_P0M];
            f[DIR_M0M]  = ftemp[DIR_M0M] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_M0M]- delf* WEIGTH[DIR_M0M];
            f[DIR_0PM]  = ftemp[DIR_0PM] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_0PM]- delf* WEIGTH[DIR_0PM];
            f[DIR_0MM]  = ftemp[DIR_0MM] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_0MM]- delf* WEIGTH[DIR_0MM];
            f[DIR_PPM] = ftemp[DIR_PPM] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_PPM]- delf* WEIGTH[DIR_PPM];
            f[DIR_MPM] = ftemp[DIR_MPM] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_MPM]- delf* WEIGTH[DIR_MPM];
            f[DIR_PMM] = ftemp[DIR_PMM] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_PMM]- delf* WEIGTH[DIR_PMM];
            f[DIR_MMM] = ftemp[DIR_MMM] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3) * f[DIR_MMM]- delf* WEIGTH[DIR_MMM];

            distributions->setDistributionInvForDirection(f[DIR_00M], x1 + DX1[DIR_00P], x2 + DX2[DIR_00P], x3 + DX3[DIR_00P], DIR_00P);
            distributions->setDistributionInvForDirection(f[DIR_P0M], x1 + DX1[DIR_M0P], x2 + DX2[DIR_M0P], x3 + DX3[DIR_M0P], DIR_M0P);
            distributions->setDistributionInvForDirection(f[DIR_M0M], x1 + DX1[DIR_P0P], x2 + DX2[DIR_P0P], x3 + DX3[DIR_P0P], DIR_P0P);
            distributions->setDistributionInvForDirection(f[DIR_0PM], x1 + DX1[DIR_0MP], x2 + DX2[DIR_0MP], x3 + DX3[DIR_0MP], DIR_0MP);
            distributions->setDistributionInvForDirection(f[DIR_0MM], x1 + DX1[DIR_0PP], x2 + DX2[DIR_0PP], x3 + DX3[DIR_0PP], DIR_0PP);
            distributions->setDistributionInvForDirection(f[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
            distributions->setDistributionInvForDirection(f[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
            distributions->setDistributionInvForDirection(f[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
            distributions->setDistributionInvForDirection(f[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
            break;
        default:
            UB_THROW(
                UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
    }
}
