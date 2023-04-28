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
//! \file MultiphasePressureBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphasePressureBCStrategy.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

MultiphasePressureBCStrategy::MultiphasePressureBCStrategy()
{
   BCStrategy::type = BCStrategy::MultiphasePressureBCStrategy;
   BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
MultiphasePressureBCStrategy::~MultiphasePressureBCStrategy()
{
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> MultiphasePressureBCStrategy::clone()
{
   SPtr<BCStrategy> bc(new MultiphasePressureBCStrategy());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void MultiphasePressureBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void MultiphasePressureBCStrategy::addDistributionsH(SPtr<DistributionArray3D> distributionsH)
{
	this->distributionsH = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void MultiphasePressureBCStrategy::addDistributionsH2(SPtr<DistributionArray3D> distributionsH)
{
    this->distributionsH2 = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void MultiphasePressureBCStrategy::applyBC()
{
   using namespace vf::lbm::dir;

   LBMReal f[D3Q27System::ENDF+1];
   LBMReal h[D3Q27System::ENDF+1];
   LBMReal h2[D3Q27System::ENDF + 1];
   LBMReal feq[D3Q27System::ENDF+1];
   //LBMReal heq[D3Q27System::ENDF+1];
   LBMReal htemp[D3Q27System::ENDF+1];
   
   distributions->getDistributionInv(f, x1, x2, x3);
   distributionsH->getDistributionInv(h, x1, x2, x3);
   if (distributionsH2)
       distributionsH2->getDistributionInv(h2, x1, x2, x3);
   LBMReal phi, vx1, vx2, vx3, p1, phiBC;
   
   D3Q27System::calcDensity(h, phi);

   calcMacrosFct(f, p1, vx1, vx2, vx3);
   p1 = 0.0;

   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;

   //flag points in direction of fluid
   if      (bcPtr->hasVelocityBoundaryFlag(DIR_P00)) { nx1 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_M00)) { nx1 += 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_0P0)) { nx2 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_0M0)) { nx2 += 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_00P)) { nx3 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_00M)) { nx3 += 1; }
   //else UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on velocity boundary..."));
   
   phiBC = bcPtr->getBoundaryPhaseField();
   LBMReal rhoBC = bcPtr->getBoundaryDensity();
   D3Q27System::calcIncompFeq(feq, rhoBC, vx1, vx2, vx3);

   D3Q27System::calcMultiphaseHeq(htemp, phiBC, vx1, vx2, vx3);

   for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++) {
       if (bcPtr->hasDensityBoundaryFlag(fdir)) {
           LBMReal ftemp = -f[D3Q27System::INVDIR[fdir]] + feq[fdir] + feq[D3Q27System::INVDIR[fdir]];
           distributions->setDistributionForDirection(ftemp, x1, x2, x3, D3Q27System::INVDIR[fdir]);

           LBMReal hReturn = -h[D3Q27System::INVDIR[fdir]] + htemp[fdir] + htemp[D3Q27System::INVDIR[fdir]];
           distributionsH->setDistributionForDirection(hReturn, x1, x2, x3, D3Q27System::INVDIR[fdir]);
       }
   }
}
