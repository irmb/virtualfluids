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
//! \file ThixotropyVelocityWithDensityBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "ThixotropyVelocityWithDensityBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BCArray3D.h"

ThixotropyVelocityWithDensityBCAlgorithm::ThixotropyVelocityWithDensityBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::ThixotropyVelocityWithDensityBCAlgorithm;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
ThixotropyVelocityWithDensityBCAlgorithm::~ThixotropyVelocityWithDensityBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> ThixotropyVelocityWithDensityBCAlgorithm::clone()
{
   SPtr<BCAlgorithm> bc(new ThixotropyVelocityWithDensityBCAlgorithm());
   dynamicPointerCast<ThixotropyVelocityWithDensityBCAlgorithm>(bc)->setLambdaBC(lambdaBC);
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyVelocityWithDensityBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyVelocityWithDensityBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
   this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyVelocityWithDensityBCAlgorithm::applyBC()
{
    using namespace vf::lbm::dir;

   //velocity bc for non reflecting pressure bc
   real f[D3Q27System::ENDF+1];
   distributions->getDistributionInv(f, x1, x2, x3);
   
   real h[D3Q27System::ENDF + 1];
   distributionsH->getDistributionInv(h, x1, x2, x3);

   real rho, vx1, vx2, vx3, drho;
   calcMacrosFct(f, drho, vx1, vx2, vx3);
   
   rho = vf::lbm::constant::c1o1+drho*compressibleFactor;
  
   ///////////////////////////////////////////////////////////////////
   // Rheology
   real lambda = D3Q27System::getDensity(h);

   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;

   //flag points in direction of fluid
   if (bcPtr->hasVelocityBoundaryFlag(DIR_P00)) { nx1 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_M00)) { nx1 += 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_0P0)) { nx2 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_0M0)) { nx2 += 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_00P)) { nx3 -= 1; }
   else if (bcPtr->hasVelocityBoundaryFlag(DIR_00M)) { nx3 += 1; }
   else	 UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on velocity boundary..."));

   for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
   {
      int nX1 = x1 + D3Q27System::DX1[fdir];
      int nX2 = x2 + D3Q27System::DX2[fdir];
      int nX3 = x3 + D3Q27System::DX3[fdir];

      int minX1 = 0;
      int minX2 = 0;
      int minX3 = 0;

      int maxX1 = (int)bcArray->getNX1();
      int maxX2 = (int)bcArray->getNX2();
      int maxX3 = (int)bcArray->getNX3();

      if (minX1 <= nX1 && maxX1 > nX1 && minX2 <= nX2 && maxX2 > nX2 && minX3 <= nX3 && maxX3 > nX3)
      {
         if (bcArray->isSolid(nX1,nX2,nX3))
         {
            const int invDir = D3Q27System::INVDIR[fdir];
            real velocity = bcPtr->getBoundaryVelocity(fdir);

            real fReturn = (f[fdir] + f[invDir] - velocity*rho) / vf::lbm::constant::c2o1 - drho*D3Q27System::WEIGTH[invDir];
            distributions->setDistributionForDirection(fReturn, nX1, nX2, nX3, invDir);
         }
      }
      
      if (bcPtr->hasVelocityBoundaryFlag(fdir))
      {
         real htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
         htemp = D3Q27System::getCompFeqForDirection(fdir, lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
         distributionsH->setDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
      }
   }
}
