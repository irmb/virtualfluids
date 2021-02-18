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
//! \file MultiphaseVelocityBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseVelocityBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

MultiphaseVelocityBCAlgorithm::MultiphaseVelocityBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::VelocityBCAlgorithm;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
MultiphaseVelocityBCAlgorithm::~MultiphaseVelocityBCAlgorithm()
{
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> MultiphaseVelocityBCAlgorithm::clone()
{
   SPtr<BCAlgorithm> bc(new MultiphaseVelocityBCAlgorithm());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributionsH)
{
	this->distributionsH = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityBCAlgorithm::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal h[D3Q27System::ENDF+1];
   LBMReal feq[D3Q27System::ENDF+1];
   LBMReal heq[D3Q27System::ENDF+1];
   LBMReal htemp[D3Q27System::ENDF+1];
   
   distributions->getDistributionInv(f, x1, x2, x3);
   distributionsH->getDistributionInv(h, x1, x2, x3);
   LBMReal phi, vx1, vx2, vx3, p1, phiBC;
   
   D3Q27System::calcDensity(h, phi);

   calcMacrosFct(f, p1, vx1, vx2, vx3);
   D3Q27System::calcMultiphaseFeqVB(feq, p1, vx1, vx2, vx3);
   D3Q27System::calcMultiphaseHeq(heq, phi, vx1, vx2, vx3);

   ///// added for phase field //////

   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;

   //flag points in direction of fluid
   if      (bcPtr->hasVelocityBoundaryFlag(D3Q27System::E)) { nx1 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::W)) { nx1 += 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::N)) { nx2 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::S)) { nx2 += 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::T)) { nx3 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(D3Q27System::B)) { nx3 += 1; }
   else UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on velocity boundary..."));
   
   phiBC = bcPtr->getBoundaryPhaseField();
   
   D3Q27System::calcMultiphaseHeq(htemp, phiBC, vx1, vx2, vx3);

   for (int fdir = D3Q27System::STARTF; fdir<=D3Q27System::ENDF; fdir++)
   {
	   if (bcPtr->hasVelocityBoundaryFlag(fdir))
	   {
		   LBMReal hReturn = htemp[fdir]+h[fdir]-heq[fdir];
		   distributionsH->setDistributionForDirection(hReturn, nx1, nx2, nx3, fdir);
	   }
   }
   
   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasVelocityBoundaryFlag(fdir))
      {
         const int invDir = D3Q27System::INVDIR[fdir];
         LBMReal q = bcPtr->getQ(invDir);// m+m q=0 stabiler
         LBMReal velocity = bcPtr->getBoundaryVelocity(invDir);
		 LBMReal fReturn = ((1.0-q)/(1.0+q))*((f[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q*(f[invDir]+f[fdir])-velocity)/(1.0+q));
         distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }

}

