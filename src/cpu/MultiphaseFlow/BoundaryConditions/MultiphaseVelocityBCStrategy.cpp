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
//! \file MultiphaseVelocityBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseVelocityBCStrategy.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

MultiphaseVelocityBCStrategy::MultiphaseVelocityBCStrategy()
{
   BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
MultiphaseVelocityBCStrategy::~MultiphaseVelocityBCStrategy()
{
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> MultiphaseVelocityBCStrategy::clone()
{
   SPtr<BCStrategy> bc(new MultiphaseVelocityBCStrategy());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityBCStrategy::addDistributionsH(SPtr<DistributionArray3D> distributionsH)
{
	this->distributionsH = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityBCStrategy::addDistributionsH2(SPtr<DistributionArray3D> distributionsH)
{
    this->distributionsH2 = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseVelocityBCStrategy::applyBC()
{
    using namespace vf::lbm::dir;

   real f[D3Q27System::ENDF+1];
   real h[D3Q27System::ENDF+1];
   real h2[D3Q27System::ENDF + 1];
   real feq[D3Q27System::ENDF+1];
   real heq[D3Q27System::ENDF+1];
   real htemp[D3Q27System::ENDF+1];
   
   distributions->getDistributionInv(f, x1, x2, x3);
   distributionsH->getDistributionInv(h, x1, x2, x3);
   if (distributionsH2)
       distributionsH2->getDistributionInv(h2, x1, x2, x3);
   real phi, vx1, vx2, vx3, p1, phiBC;
   
   D3Q27System::calcDensity(h, phi);

   calcMacrosFct(f, p1, vx1, vx2, vx3);
   vx1=bcPtr->getBoundaryVelocityX1();
   vx2 = bcPtr->getBoundaryVelocityX2();
   vx3 = bcPtr->getBoundaryVelocityX3();
   p1 = vf::basics::constant::c0o1;
   D3Q27System::calcMultiphaseFeqVB(feq, p1, vx1, vx2, vx3);
   D3Q27System::calcMultiphaseHeq(heq, phi, vx1, vx2, vx3);

   ///// added for phase field //////

   phiBC = bcPtr->getBoundaryPhaseField();
   
   D3Q27System::calcMultiphaseHeq(htemp, phiBC, vx1, vx2, vx3);
   //D3Q27System::calcMultiphaseHeq(htemp, phiBC, bcPtr->getBoundaryVelocityX1(), bcPtr->getBoundaryVelocityX2(), bcPtr->getBoundaryVelocityX2());//30.03.2021 EQ phase field BC!
   //for (int fdir = D3Q27System::STARTF; fdir<=D3Q27System::ENDF; fdir++)
   //{
	  // if (bcPtr->hasVelocityBoundaryFlag(fdir))
	  // {
		 //  LBMReal hReturn = htemp[fdir]+h[fdir]-heq[fdir];
   //        //17.03.2021 Let us just set the plain eq
   //        //LBMReal hReturn = htemp[fdir];
		 //  distributionsH->setDistributionForDirection(hReturn, nx1, nx2, nx3, fdir);
   //      //  if (distributionsH2)
   //      //      distributionsH2->setDistributionForDirection(0, nx1, nx2, nx3, fdir);
	  // }
   //}
   
   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasVelocityBoundaryFlag(fdir))
      {
         const int invDir = D3Q27System::INVDIR[fdir];
         //LBMReal q = bcPtr->getQ(invDir);// m+m q=0 stabiler
         real velocity = bcPtr->getBoundaryVelocity(invDir);
		 //16.03.2021 quick fix for velocity BC
         real fReturn = f[invDir] - velocity;
         //LBMReal fReturn = ((1.0-q)/(1.0+q))*((f[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q*(f[invDir]+f[fdir])-velocity)/(1.0+q));
        // distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);//no delay BB
         distributions->setDistributionForDirection(fReturn, x1, x2, x3, invDir);//delay BB  

         real hReturn = htemp[invDir] + h[invDir] - heq[invDir] - velocity*phi;
         distributionsH->setDistributionForDirection(hReturn, x1, x2, x3, invDir);//delay BB  
         if (distributionsH2) {
             fReturn = h2[invDir] ;
            // distributionsH2->setDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
             distributionsH2->setDistributionForDirection(fReturn, x1, x2, x3, invDir);//delay BB 
         }

      }
   }

}

