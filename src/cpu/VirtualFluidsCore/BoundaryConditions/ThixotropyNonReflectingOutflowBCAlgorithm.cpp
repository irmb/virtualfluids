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
//! \file ThixotropyNonReflectingOutflowBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "ThixotropyNonReflectingOutflowBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

ThixotropyNonReflectingOutflowBCAlgorithm::ThixotropyNonReflectingOutflowBCAlgorithm()
{
	BCAlgorithm::type = BCAlgorithm::ThixotropyNonReflectingOutflowBCAlgorithm;
	BCAlgorithm::preCollision = true;
	BCAlgorithm::thixotropy = true;
}
//////////////////////////////////////////////////////////////////////////
ThixotropyNonReflectingOutflowBCAlgorithm::~ThixotropyNonReflectingOutflowBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> ThixotropyNonReflectingOutflowBCAlgorithm::clone()
{
	SPtr<BCAlgorithm> bc(new ThixotropyNonReflectingOutflowBCAlgorithm());
	return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyNonReflectingOutflowBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
	this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
//void ThixotropyNonReflectingOutflowBCAlgorithm::addDistributionsF(DistributionArray3DPtr distributions)
//{
//	this->distributionsf = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void ThixotropyNonReflectingOutflowBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
	this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyNonReflectingOutflowBCAlgorithm::applyBC()
{
   using namespace D3Q27System;
   LBMReal f[ENDF + 1];
   LBMReal ftemp[ENDF + 1];

   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;
   int direction = -1;

   //flag points in direction of fluid
   if (bcPtr->hasDensityBoundaryFlag(E)) { nx1 += 1; direction = E; }
   else if (bcPtr->hasDensityBoundaryFlag(W)) { nx1 -= 1; direction = W; }
   else if (bcPtr->hasDensityBoundaryFlag(N)) { nx2 += 1; direction = N; }
   else if (bcPtr->hasDensityBoundaryFlag(S)) { nx2 -= 1; direction = S; }
   else if (bcPtr->hasDensityBoundaryFlag(T)) { nx3 += 1; direction = T; }
   else if (bcPtr->hasDensityBoundaryFlag(B)) { nx3 -= 1; direction = B; }
   else UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

   distributions->getDistribution(f, x1, x2, x3);
   distributions->getDistribution(ftemp, nx1, nx2, nx3);

   LBMReal rho, vx1, vx2, vx3;
   calcMacrosFct(f, rho, vx1, vx2, vx3);

   switch (direction)
   {
   case E:
      f[DIR_P00] = ftemp[DIR_P00] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[DIR_P00];
      f[NE] = ftemp[NE] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[NE];
      f[DIR_PM0] = ftemp[DIR_PM0] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[DIR_PM0];
      f[DIR_P0P] = ftemp[DIR_P0P] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[DIR_P0P];
      f[DIR_P0M] = ftemp[DIR_P0M] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[DIR_P0M];
      f[DIR_PPP] = ftemp[DIR_PPP] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[DIR_PPP];
      f[DIR_PMP] = ftemp[DIR_PMP] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[DIR_PMP];
      f[DIR_PPM] = ftemp[DIR_PPM] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[DIR_PPM];
      f[DIR_PMM] = ftemp[DIR_PMM] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * f[DIR_PMM];

      distributions->setDistributionInvForDirection(f[DIR_P00], x1 + DX1[W], x2 + DX2[W], x3 + DX3[W], W);
      distributions->setDistributionInvForDirection(f[NE], x1 + DX1[DIR_MM0], x2 + DX2[DIR_MM0], x3 + DX3[DIR_MM0], DIR_MM0);
      distributions->setDistributionInvForDirection(f[DIR_PM0], x1 + DX1[DIR_MP0], x2 + DX2[DIR_MP0], x3 + DX3[DIR_MP0], DIR_MP0);
      distributions->setDistributionInvForDirection(f[DIR_P0P], x1 + DX1[DIR_M0M], x2 + DX2[DIR_M0M], x3 + DX3[DIR_M0M], DIR_M0M);
      distributions->setDistributionInvForDirection(f[DIR_P0M], x1 + DX1[DIR_M0P], x2 + DX2[DIR_M0P], x3 + DX3[DIR_M0P], DIR_M0P);
      distributions->setDistributionInvForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
      distributions->setDistributionInvForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
      distributions->setDistributionInvForDirection(f[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
      distributions->setDistributionInvForDirection(f[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
      break;
   case W:
      f[W] = ftemp[W] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[W];
      f[DIR_MP0] = ftemp[DIR_MP0] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[DIR_MP0];
      f[DIR_MM0] = ftemp[DIR_MM0] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[DIR_MM0];
      f[DIR_M0P] = ftemp[DIR_M0P] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[DIR_M0P];
      f[DIR_M0M] = ftemp[DIR_M0M] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[DIR_M0M];
      f[DIR_MPP] = ftemp[DIR_MPP] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[DIR_MPP];
      f[DIR_MMP] = ftemp[DIR_MMP] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[DIR_MMP];
      f[DIR_MPM] = ftemp[DIR_MPM] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[DIR_MPM];
      f[DIR_MMM] = ftemp[DIR_MMM] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * f[DIR_MMM];

      distributions->setDistributionInvForDirection(f[W], x1 + DX1[DIR_P00], x2 + DX2[DIR_P00], x3 + DX3[DIR_P00], E);
      distributions->setDistributionInvForDirection(f[DIR_MP0], x1 + DX1[DIR_PM0], x2 + DX2[DIR_PM0], x3 + DX3[DIR_PM0], DIR_PM0);
      distributions->setDistributionInvForDirection(f[DIR_MM0], x1 + DX1[NE], x2 + DX2[NE], x3 + DX3[NE], NE);
      distributions->setDistributionInvForDirection(f[DIR_M0P], x1 + DX1[DIR_P0M], x2 + DX2[DIR_P0M], x3 + DX3[DIR_P0M], DIR_P0M);
      distributions->setDistributionInvForDirection(f[DIR_M0M], x1 + DX1[DIR_P0P], x2 + DX2[DIR_P0P], x3 + DX3[DIR_P0P], DIR_P0P);
      distributions->setDistributionInvForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
      distributions->setDistributionInvForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
      distributions->setDistributionInvForDirection(f[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
      distributions->setDistributionInvForDirection(f[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
      break;
   case N:
      f[N] = ftemp[N] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[N];
      f[NE] = ftemp[NE] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[NE];
      f[DIR_MP0] = ftemp[DIR_MP0] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[DIR_MP0];
      f[DIR_0PP] = ftemp[DIR_0PP] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[DIR_0PP];
      f[DIR_0PM] = ftemp[DIR_0PM] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[DIR_0PM];
      f[DIR_PPP] = ftemp[DIR_PPP] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[DIR_PPP];
      f[DIR_MPP] = ftemp[DIR_MPP] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[DIR_MPP];
      f[DIR_PPM] = ftemp[DIR_PPM] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[DIR_PPM];
      f[DIR_MPM] = ftemp[DIR_MPM] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * f[DIR_MPM];

      distributions->setDistributionInvForDirection(f[N], x1 + DX1[S], x2 + DX2[S], x3 + DX3[S], S);
      distributions->setDistributionInvForDirection(f[NE], x1 + DX1[DIR_MM0], x2 + DX2[DIR_MM0], x3 + DX3[DIR_MM0], DIR_MM0);
      distributions->setDistributionInvForDirection(f[DIR_MP0], x1 + DX1[DIR_PM0], x2 + DX2[DIR_PM0], x3 + DX3[DIR_PM0], DIR_PM0);
      distributions->setDistributionInvForDirection(f[DIR_0PP], x1 + DX1[DIR_0MM], x2 + DX2[DIR_0MM], x3 + DX3[DIR_0MM], DIR_0MM);
      distributions->setDistributionInvForDirection(f[DIR_0PM], x1 + DX1[DIR_0MP], x2 + DX2[DIR_0MP], x3 + DX3[DIR_0MP], DIR_0MP);
      distributions->setDistributionInvForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
      distributions->setDistributionInvForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
      distributions->setDistributionInvForDirection(f[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
      distributions->setDistributionInvForDirection(f[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
      break;
   case S:
      f[S] = ftemp[S] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[S];
      f[DIR_PM0] = ftemp[DIR_PM0] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[DIR_PM0];
      f[DIR_MM0] = ftemp[DIR_MM0] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[DIR_MM0];
      f[DIR_0MP] = ftemp[DIR_0MP] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[DIR_0MP];
      f[DIR_0MM] = ftemp[DIR_0MM] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[DIR_0MM];
      f[DIR_PMP] = ftemp[DIR_PMP] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[DIR_PMP];
      f[DIR_MMP] = ftemp[DIR_MMP] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[DIR_MMP];
      f[DIR_PMM] = ftemp[DIR_PMM] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[DIR_PMM];
      f[DIR_MMM] = ftemp[DIR_MMM] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * f[DIR_MMM];

      distributions->setDistributionInvForDirection(f[S], x1 + DX1[N], x2 + DX2[N], x3 + DX3[N], N);
      distributions->setDistributionInvForDirection(f[DIR_PM0], x1 + DX1[DIR_MP0], x2 + DX2[DIR_MP0], x3 + DX3[DIR_MP0], DIR_MP0);
      distributions->setDistributionInvForDirection(f[DIR_MM0], x1 + DX1[NE], x2 + DX2[NE], x3 + DX3[NE], NE);
      distributions->setDistributionInvForDirection(f[DIR_0MP], x1 + DX1[DIR_0PM], x2 + DX2[DIR_0PM], x3 + DX3[DIR_0PM], DIR_0PM);
      distributions->setDistributionInvForDirection(f[DIR_0MM], x1 + DX1[DIR_0PP], x2 + DX2[DIR_0PP], x3 + DX3[DIR_0PP], DIR_0PP);
      distributions->setDistributionInvForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
      distributions->setDistributionInvForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
      distributions->setDistributionInvForDirection(f[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
      distributions->setDistributionInvForDirection(f[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
      break;
   case T:
      f[T] = ftemp[T] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[T];
      f[DIR_P0P] = ftemp[DIR_P0P] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[DIR_P0P];
      f[DIR_M0P] = ftemp[DIR_M0P] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[DIR_M0P];
      f[DIR_0PP] = ftemp[DIR_0PP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[DIR_0PP];
      f[DIR_0MP] = ftemp[DIR_0MP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[DIR_0MP];
      f[DIR_PPP] = ftemp[DIR_PPP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[DIR_PPP];
      f[DIR_MPP] = ftemp[DIR_MPP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[DIR_MPP];
      f[DIR_PMP] = ftemp[DIR_PMP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[DIR_PMP];
      f[DIR_MMP] = ftemp[DIR_MMP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * f[DIR_MMP];

      distributions->setDistributionInvForDirection(f[T], x1 + DX1[B], x2 + DX2[B], x3 + DX3[B], B);
      distributions->setDistributionInvForDirection(f[DIR_P0P], x1 + DX1[DIR_M0M], x2 + DX2[DIR_M0M], x3 + DX3[DIR_M0M], DIR_M0M);
      distributions->setDistributionInvForDirection(f[DIR_M0P], x1 + DX1[DIR_P0M], x2 + DX2[DIR_P0M], x3 + DX3[DIR_P0M], DIR_P0M);
      distributions->setDistributionInvForDirection(f[DIR_0PP], x1 + DX1[DIR_0MM], x2 + DX2[DIR_0MM], x3 + DX3[DIR_0MM], DIR_0MM);
      distributions->setDistributionInvForDirection(f[DIR_0MP], x1 + DX1[DIR_0PM], x2 + DX2[DIR_0PM], x3 + DX3[DIR_0PM], DIR_0PM);
      distributions->setDistributionInvForDirection(f[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
      distributions->setDistributionInvForDirection(f[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
      distributions->setDistributionInvForDirection(f[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
      distributions->setDistributionInvForDirection(f[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
      break;
   case B:
      f[B] = ftemp[B] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[B];
      f[DIR_P0M] = ftemp[DIR_P0M] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[DIR_P0M];
      f[DIR_M0M] = ftemp[DIR_M0M] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[DIR_M0M];
      f[DIR_0PM] = ftemp[DIR_0PM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[DIR_0PM];
      f[DIR_0MM] = ftemp[DIR_0MM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[DIR_0MM];
      f[DIR_PPM] = ftemp[DIR_PPM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[DIR_PPM];
      f[DIR_MPM] = ftemp[DIR_MPM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[DIR_MPM];
      f[DIR_PMM] = ftemp[DIR_PMM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[DIR_PMM];
      f[DIR_MMM] = ftemp[DIR_MMM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * f[DIR_MMM];

      distributions->setDistributionInvForDirection(f[B], x1 + DX1[T], x2 + DX2[T], x3 + DX3[T], T);
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
      UB_THROW(UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
   }
   LBMReal h[D3Q27System::ENDF + 1];
   LBMReal htemp[ENDF + 1];

   distributionsH->getDistribution(h, x1, x2, x3);
   distributionsH->getDistribution(htemp, nx1, nx2, nx3);

   vx1 = 0.0;
   vx2 = 0.0;
   vx3 = 0.0;

   //LBMReal rho, vx1, vx2, vx3;
   //calcMacrosFct(f, rho, vx1, vx2, vx3);

   switch (direction)
   {
   case E:
      h[DIR_P00]  = htemp[DIR_P00] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[DIR_P00];
      h[NE] = htemp[NE] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[NE];
      h[DIR_PM0] = htemp[DIR_PM0] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[DIR_PM0];
      h[DIR_P0P] = htemp[DIR_P0P] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[DIR_P0P];
      h[DIR_P0M] = htemp[DIR_P0M] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[DIR_P0M];
      h[DIR_PPP] = htemp[DIR_PPP] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[DIR_PPP];
      h[DIR_PMP] = htemp[DIR_PMP] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[DIR_PMP];
      h[DIR_PPM] = htemp[DIR_PPM] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[DIR_PPM];
      h[DIR_PMM] = htemp[DIR_PMM] * (UbMath::one_over_sqrt3 + vx1) + (1.0 - UbMath::one_over_sqrt3 - vx1) * h[DIR_PMM];

      distributionsH->setDistributionInvForDirection(h[DIR_P00], x1 + DX1[W], x2 + DX2[W], x3 + DX3[W], W);
      distributionsH->setDistributionInvForDirection(h[NE], x1 + DX1[DIR_MM0], x2 + DX2[DIR_MM0], x3 + DX3[DIR_MM0], DIR_MM0);
      distributionsH->setDistributionInvForDirection(h[DIR_PM0], x1 + DX1[DIR_MP0], x2 + DX2[DIR_MP0], x3 + DX3[DIR_MP0], DIR_MP0);
      distributionsH->setDistributionInvForDirection(h[DIR_P0P], x1 + DX1[DIR_M0M], x2 + DX2[DIR_M0M], x3 + DX3[DIR_M0M], DIR_M0M);
      distributionsH->setDistributionInvForDirection(h[DIR_P0M], x1 + DX1[DIR_M0P], x2 + DX2[DIR_M0P], x3 + DX3[DIR_M0P], DIR_M0P);
      distributionsH->setDistributionInvForDirection(h[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
      distributionsH->setDistributionInvForDirection(h[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
      distributionsH->setDistributionInvForDirection(h[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
      distributionsH->setDistributionInvForDirection(h[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
      break;
   case W:
      h[W] = htemp[W] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[W];
      h[DIR_MP0] = htemp[DIR_MP0] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[DIR_MP0];
      h[DIR_MM0] = htemp[DIR_MM0] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[DIR_MM0];
      h[DIR_M0P] = htemp[DIR_M0P] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[DIR_M0P];
      h[DIR_M0M] = htemp[DIR_M0M] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[DIR_M0M];
      h[DIR_MPP] = htemp[DIR_MPP] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[DIR_MPP];
      h[DIR_MMP] = htemp[DIR_MMP] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[DIR_MMP];
      h[DIR_MPM] = htemp[DIR_MPM] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[DIR_MPM];
      h[DIR_MMM] = htemp[DIR_MMM] * (UbMath::one_over_sqrt3 - vx1) + (1.0 - UbMath::one_over_sqrt3 + vx1) * h[DIR_MMM];

      distributionsH->setDistributionInvForDirection(h[W], x1 + DX1[DIR_P00], x2 + DX2[DIR_P00], x3 + DX3[DIR_P00], E);
      distributionsH->setDistributionInvForDirection(h[DIR_MP0], x1 + DX1[DIR_PM0], x2 + DX2[DIR_PM0], x3 + DX3[DIR_PM0], DIR_PM0);
      distributionsH->setDistributionInvForDirection(h[DIR_MM0], x1 + DX1[NE], x2 + DX2[NE], x3 + DX3[NE], NE);
      distributionsH->setDistributionInvForDirection(h[DIR_M0P], x1 + DX1[DIR_P0M], x2 + DX2[DIR_P0M], x3 + DX3[DIR_P0M], DIR_P0M);
      distributionsH->setDistributionInvForDirection(h[DIR_M0M], x1 + DX1[DIR_P0P], x2 + DX2[DIR_P0P], x3 + DX3[DIR_P0P], DIR_P0P);
      distributionsH->setDistributionInvForDirection(h[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
      distributionsH->setDistributionInvForDirection(h[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
      distributionsH->setDistributionInvForDirection(h[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
      distributionsH->setDistributionInvForDirection(h[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
      break;
   case N:
      h[N] = htemp[N] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[N];
      h[NE] = htemp[NE] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[NE];
      h[DIR_MP0] = htemp[DIR_MP0] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[DIR_MP0];
      h[DIR_0PP] = htemp[DIR_0PP] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[DIR_0PP];
      h[DIR_0PM] = htemp[DIR_0PM] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[DIR_0PM];
      h[DIR_PPP] = htemp[DIR_PPP] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[DIR_PPP];
      h[DIR_MPP] = htemp[DIR_MPP] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[DIR_MPP];
      h[DIR_PPM] = htemp[DIR_PPM] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[DIR_PPM];
      h[DIR_MPM] = htemp[DIR_MPM] * (UbMath::one_over_sqrt3 + vx2) + (1.0 - UbMath::one_over_sqrt3 - vx2) * h[DIR_MPM];

      distributionsH->setDistributionInvForDirection(h[N], x1 + DX1[S], x2 + DX2[S], x3 + DX3[S], S);
      distributionsH->setDistributionInvForDirection(h[NE], x1 + DX1[DIR_MM0], x2 + DX2[DIR_MM0], x3 + DX3[DIR_MM0], DIR_MM0);
      distributionsH->setDistributionInvForDirection(h[DIR_MP0], x1 + DX1[DIR_PM0], x2 + DX2[DIR_PM0], x3 + DX3[DIR_PM0], DIR_PM0);
      distributionsH->setDistributionInvForDirection(h[DIR_0PP], x1 + DX1[DIR_0MM], x2 + DX2[DIR_0MM], x3 + DX3[DIR_0MM], DIR_0MM);
      distributionsH->setDistributionInvForDirection(h[DIR_0PM], x1 + DX1[DIR_0MP], x2 + DX2[DIR_0MP], x3 + DX3[DIR_0MP], DIR_0MP);
      distributionsH->setDistributionInvForDirection(h[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
      distributionsH->setDistributionInvForDirection(h[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
      distributionsH->setDistributionInvForDirection(h[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
      distributionsH->setDistributionInvForDirection(h[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
      break;
   case S:
      h[S] = htemp[S] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[S];
      h[DIR_PM0] = htemp[DIR_PM0] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[DIR_PM0];
      h[DIR_MM0] = htemp[DIR_MM0] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[DIR_MM0];
      h[DIR_0MP] = htemp[DIR_0MP] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[DIR_0MP];
      h[DIR_0MM] = htemp[DIR_0MM] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[DIR_0MM];
      h[DIR_PMP] = htemp[DIR_PMP] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[DIR_PMP];
      h[DIR_MMP] = htemp[DIR_MMP] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[DIR_MMP];
      h[DIR_PMM] = htemp[DIR_PMM] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[DIR_PMM];
      h[DIR_MMM] = htemp[DIR_MMM] * (UbMath::one_over_sqrt3 - vx2) + (1.0 - UbMath::one_over_sqrt3 + vx2) * h[DIR_MMM];

      distributionsH->setDistributionInvForDirection(h[S], x1 + DX1[N], x2 + DX2[N], x3 + DX3[N], N);
      distributionsH->setDistributionInvForDirection(h[DIR_PM0], x1 + DX1[DIR_MP0], x2 + DX2[DIR_MP0], x3 + DX3[DIR_MP0], DIR_MP0);
      distributionsH->setDistributionInvForDirection(h[DIR_MM0], x1 + DX1[NE], x2 + DX2[NE], x3 + DX3[NE], NE);
      distributionsH->setDistributionInvForDirection(h[DIR_0MP], x1 + DX1[DIR_0PM], x2 + DX2[DIR_0PM], x3 + DX3[DIR_0PM], DIR_0PM);
      distributionsH->setDistributionInvForDirection(h[DIR_0MM], x1 + DX1[DIR_0PP], x2 + DX2[DIR_0PP], x3 + DX3[DIR_0PP], DIR_0PP);
      distributionsH->setDistributionInvForDirection(h[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
      distributionsH->setDistributionInvForDirection(h[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
      distributionsH->setDistributionInvForDirection(h[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
      distributionsH->setDistributionInvForDirection(h[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
      break;
   case T:
      h[T] = htemp[T] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[T];
      h[DIR_P0P] = htemp[DIR_P0P] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[DIR_P0P];
      h[DIR_M0P] = htemp[DIR_M0P] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[DIR_M0P];
      h[DIR_0PP] = htemp[DIR_0PP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[DIR_0PP];
      h[DIR_0MP] = htemp[DIR_0MP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[DIR_0MP];
      h[DIR_PPP] = htemp[DIR_PPP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[DIR_PPP];
      h[DIR_MPP] = htemp[DIR_MPP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[DIR_MPP];
      h[DIR_PMP] = htemp[DIR_PMP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[DIR_PMP];
      h[DIR_MMP] = htemp[DIR_MMP] * (UbMath::one_over_sqrt3 + vx3) + (1.0 - UbMath::one_over_sqrt3 - vx3) * h[DIR_MMP];

      distributionsH->setDistributionInvForDirection(h[T], x1 + DX1[B], x2 + DX2[B], x3 + DX3[B], B);
      distributionsH->setDistributionInvForDirection(h[DIR_P0P], x1 + DX1[DIR_M0M], x2 + DX2[DIR_M0M], x3 + DX3[DIR_M0M], DIR_M0M);
      distributionsH->setDistributionInvForDirection(h[DIR_M0P], x1 + DX1[DIR_P0M], x2 + DX2[DIR_P0M], x3 + DX3[DIR_P0M], DIR_P0M);
      distributionsH->setDistributionInvForDirection(h[DIR_0PP], x1 + DX1[DIR_0MM], x2 + DX2[DIR_0MM], x3 + DX3[DIR_0MM], DIR_0MM);
      distributionsH->setDistributionInvForDirection(h[DIR_0MP], x1 + DX1[DIR_0PM], x2 + DX2[DIR_0PM], x3 + DX3[DIR_0PM], DIR_0PM);
      distributionsH->setDistributionInvForDirection(h[DIR_PPP], x1 + DX1[DIR_MMM], x2 + DX2[DIR_MMM], x3 + DX3[DIR_MMM], DIR_MMM);
      distributionsH->setDistributionInvForDirection(h[DIR_MPP], x1 + DX1[DIR_PMM], x2 + DX2[DIR_PMM], x3 + DX3[DIR_PMM], DIR_PMM);
      distributionsH->setDistributionInvForDirection(h[DIR_PMP], x1 + DX1[DIR_MPM], x2 + DX2[DIR_MPM], x3 + DX3[DIR_MPM], DIR_MPM);
      distributionsH->setDistributionInvForDirection(h[DIR_MMP], x1 + DX1[DIR_PPM], x2 + DX2[DIR_PPM], x3 + DX3[DIR_PPM], DIR_PPM);
      break;
   case B:
      h[B] = htemp[B] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[B];
      h[DIR_P0M] = htemp[DIR_P0M] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[DIR_P0M];
      h[DIR_M0M] = htemp[DIR_M0M] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[DIR_M0M];
      h[DIR_0PM] = htemp[DIR_0PM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[DIR_0PM];
      h[DIR_0MM] = htemp[DIR_0MM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[DIR_0MM];
      h[DIR_PPM] = htemp[DIR_PPM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[DIR_PPM];
      h[DIR_MPM] = htemp[DIR_MPM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[DIR_MPM];
      h[DIR_PMM] = htemp[DIR_PMM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[DIR_PMM];
      h[DIR_MMM] = htemp[DIR_MMM] * (UbMath::one_over_sqrt3 - vx3) + (1.0 - UbMath::one_over_sqrt3 + vx3) * h[DIR_MMM];

      distributionsH->setDistributionInvForDirection(h[B], x1 + DX1[T], x2 + DX2[T], x3 + DX3[T], T);
      distributionsH->setDistributionInvForDirection(h[DIR_P0M], x1 + DX1[DIR_M0P], x2 + DX2[DIR_M0P], x3 + DX3[DIR_M0P], DIR_M0P);
      distributionsH->setDistributionInvForDirection(h[DIR_M0M], x1 + DX1[DIR_P0P], x2 + DX2[DIR_P0P], x3 + DX3[DIR_P0P], DIR_P0P);
      distributionsH->setDistributionInvForDirection(h[DIR_0PM], x1 + DX1[DIR_0MP], x2 + DX2[DIR_0MP], x3 + DX3[DIR_0MP], DIR_0MP);
      distributionsH->setDistributionInvForDirection(h[DIR_0MM], x1 + DX1[DIR_0PP], x2 + DX2[DIR_0PP], x3 + DX3[DIR_0PP], DIR_0PP);
      distributionsH->setDistributionInvForDirection(h[DIR_PPM], x1 + DX1[DIR_MMP], x2 + DX2[DIR_MMP], x3 + DX3[DIR_MMP], DIR_MMP);
      distributionsH->setDistributionInvForDirection(h[DIR_MPM], x1 + DX1[DIR_PMP], x2 + DX2[DIR_PMP], x3 + DX3[DIR_PMP], DIR_PMP);
      distributionsH->setDistributionInvForDirection(h[DIR_PMM], x1 + DX1[DIR_MPP], x2 + DX2[DIR_MPP], x3 + DX3[DIR_MPP], DIR_MPP);
      distributionsH->setDistributionInvForDirection(h[DIR_MMM], x1 + DX1[DIR_PPP], x2 + DX2[DIR_PPP], x3 + DX3[DIR_PPP], DIR_PPP);
      break;
   default:
      UB_THROW(UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
   }
}
