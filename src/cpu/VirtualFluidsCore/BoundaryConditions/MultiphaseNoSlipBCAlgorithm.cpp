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
//! \file MultiphaseNoSlipBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseNoSlipBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

MultiphaseNoSlipBCAlgorithm::MultiphaseNoSlipBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::MultiphaseNoSlipBCAlgorithm;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
MultiphaseNoSlipBCAlgorithm::~MultiphaseNoSlipBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> MultiphaseNoSlipBCAlgorithm::clone()
{
   SPtr<BCAlgorithm> bc(new MultiphaseNoSlipBCAlgorithm());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNoSlipBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNoSlipBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributionsH)
{
	this->distributionsH = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNoSlipBCAlgorithm::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal h[D3Q27System::ENDF+1];
   LBMReal h2[D3Q27System::ENDF + 1];
   //LBMReal feq[D3Q27System::ENDF+1];
   //LBMReal heq[D3Q27System::ENDF+1];
   distributions ->getDistributionInv(f, x1, x2, x3);
   if (distributionsH2)
       distributionsH2->getDistributionInv(h2, x1, x2, x3);
   distributionsH->getDistributionInv(h, x1, x2, x3);
  // LBMReal phi, vx1, vx2, vx3, p1;
   
 //  D3Q27System::calcDensity(h, phi);
   
 //  calcMacrosFct(f, p1, vx1, vx2, vx3);
 //  D3Q27System::calcMultiphaseFeqVB(feq, p1, vx1, vx2, vx3);
 //  D3Q27System::calcMultiphaseHeq(heq, phi, vx1, vx2, vx3); 

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasNoSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fdir];
		 LBMReal fReturn = f[invDir];
         //distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
         distributions->setDistributionForDirection(fReturn, x1, x2, x3, invDir);//delay BB 
         LBMReal hReturn = h[invDir];
		// distributionsH->setDistributionForDirection(hReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
         distributionsH->setDistributionForDirection(hReturn, x1, x2, x3, invDir);//delay BB  
         if (distributionsH2)
         {
             LBMReal h2Return = h2[invDir];
             distributionsH2->setDistributionForDirection(h2Return, x1, x2, x3, invDir);//delay BB
            // distributionsH2->setDistributionForDirection(h2Return, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);

         }
      }
   }
}
