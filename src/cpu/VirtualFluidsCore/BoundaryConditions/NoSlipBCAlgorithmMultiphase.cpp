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
//! \file NoSlipBCAlgorithmMultiphase.cpp
//! \ingroup BoundarConditions
//! \author Hesameddin Safari
//=======================================================================================

#include "NoSlipBCAlgorithmMultiphase.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

NoSlipBCAlgorithmMultiphase::NoSlipBCAlgorithmMultiphase()
{
   BCAlgorithm::type = BCAlgorithm::NoSlipBCAlgorithmMultiphase;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
NoSlipBCAlgorithmMultiphase::~NoSlipBCAlgorithmMultiphase()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> NoSlipBCAlgorithmMultiphase::clone()
{
   SPtr<BCAlgorithm> bc(new NoSlipBCAlgorithmMultiphase());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void NoSlipBCAlgorithmMultiphase::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void NoSlipBCAlgorithmMultiphase::addDistributionsH(SPtr<DistributionArray3D> distributionsH)
{
	this->distributionsH = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void NoSlipBCAlgorithmMultiphase::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal h[D3Q27System::ENDF+1];
   LBMReal feq[D3Q27System::ENDF+1];
   LBMReal heq[D3Q27System::ENDF+1];
   distributions ->getDistributionInv(f, x1, x2, x3);
   distributionsH->getDistributionInv(h, x1, x2, x3);
   LBMReal phi, vx1, vx2, vx3, p1;
   
   D3Q27System::calcDensity(h, phi);
   
   //LBMReal collFactorM = phi*collFactorL + (1-phi)*collFactorG;
   //LBMReal collFactorM = collFactorL + (collFactorL - collFactorG)*(phi - phiH)/(phiH - phiL);
   
   //rho = phi + (1.0 - phi)*1.0/densityRatio;
   //LBMReal rhoH = 1.0;
   //LBMReal rhoL = 1.0/densityRatio;
   //rho = rhoH + (rhoH - rhoL)*(phi - phiH)/(phiH - phiL);
   
   calcMacrosFct(f, p1, vx1, vx2, vx3);
   /*vx1/=(rho*c1o3);
   vx2/=(rho*c1o3);
   vx3/=(rho*c1o3);*/ 
   
   //calcFeqFct(feq, rho, vx1, vx2, vx3);
   //D3Q27System::calcMultiphaseFeq(feq, rho, p1, vx1, vx2, vx3);
   D3Q27System::calcMultiphaseFeqVB(feq, p1, vx1, vx2, vx3);
   D3Q27System::calcMultiphaseHeq(heq, phi, vx1, vx2, vx3); 
   //LBMReal collFactorM1 = 0.9;
   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasNoSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fdir];
         //LBMReal q = bcPtr->getQ(invDir);
         //LBMReal fReturn = ((1.0-q)/(1.0+q))*((f[invDir]-feq[invDir])/(1.0-collFactorM)+feq[invDir])+((q/(1.0+q))*(f[invDir]+f[fdir]));
		 LBMReal fReturn = f[invDir];
         distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);

		 //LBMReal hReturn = ((1.0-q)/(1.0+q))*((h[invDir]-heq[invDir])/(1.0-collFactorPh)+heq[invDir])+((q/(1.0+q))*(h[invDir]+h[fdir]));
		 LBMReal hReturn = h[invDir];
		 distributionsH->setDistributionForDirection(hReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }
}
