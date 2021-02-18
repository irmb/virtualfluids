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
//! \file MultiphaseNonReflectingOutflowBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseNonReflectingOutflowBCAlgorithm.h"
#include "D3Q27System.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

MultiphaseNonReflectingOutflowBCAlgorithm::MultiphaseNonReflectingOutflowBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::NonReflectingOutflowBCAlgorithm;
   BCAlgorithm::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
MultiphaseNonReflectingOutflowBCAlgorithm::~MultiphaseNonReflectingOutflowBCAlgorithm()
{
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> MultiphaseNonReflectingOutflowBCAlgorithm::clone()
{
    SPtr<BCAlgorithm> bc(new MultiphaseNonReflectingOutflowBCAlgorithm());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNonReflectingOutflowBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNonReflectingOutflowBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributionsH)
{
	this->distributionsH = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNonReflectingOutflowBCAlgorithm::applyBC()
{
   using namespace D3Q27System;
   using namespace UbMath;
   LBMReal f[ENDF+1];
   LBMReal ftemp[ENDF+1];
   LBMReal h[D3Q27System::ENDF+1];
   LBMReal htemp[ENDF+1];

   int nx1 = x1;
   int nx2 = x2;
   int nx3 = x3;
   int direction = -1;

   //flag points in direction of fluid
   if      (bcPtr->hasDensityBoundaryFlag(E)) { nx1 += 1; direction = E; }
   else if (bcPtr->hasDensityBoundaryFlag(W)) { nx1 -= 1; direction = W; }
   else if (bcPtr->hasDensityBoundaryFlag(N)) { nx2 += 1; direction = N; }
   else if (bcPtr->hasDensityBoundaryFlag(S)) { nx2 -= 1; direction = S; }
   else if (bcPtr->hasDensityBoundaryFlag(T)) { nx3 += 1; direction = T; }
   else if (bcPtr->hasDensityBoundaryFlag(B)) { nx3 -= 1; direction = B; }
   else UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

   distributions->getDistribution(f, x1, x2, x3);
   distributions->getDistribution(ftemp, nx1, nx2, nx3);
   distributionsH->getDistribution(h, x1, x2, x3);
   distributionsH->getDistribution(htemp, nx1, nx2, nx3);

   LBMReal phi, p1, vx1, vx2, vx3;
   
   D3Q27System::calcDensity(h, phi);
   
   calcMacrosFct(f, p1, vx1, vx2, vx3);

   switch (direction)
   {
   case E:
      f[E]   = ftemp[E]   * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[E]   ;
      f[NE]  = ftemp[NE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[NE]  ;
      f[SE]  = ftemp[SE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[SE]  ;
      f[TE]  = ftemp[TE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[TE]  ;
      f[BE]  = ftemp[BE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[BE]  ;
      f[TNE] = ftemp[TNE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[TNE] ;
      f[TSE] = ftemp[TSE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[TSE] ;
      f[BNE] = ftemp[BNE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[BNE] ;
      f[BSE] = ftemp[BSE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*f[BSE] ;

      distributions->setDistributionInvForDirection(f[E],   x1+DX1[W],   x2+DX2[W],   x3+DX3[W],   W);
      distributions->setDistributionInvForDirection(f[NE],  x1+DX1[SW],  x2+DX2[SW],  x3+DX3[SW],  SW);
      distributions->setDistributionInvForDirection(f[SE],  x1+DX1[NW],  x2+DX2[NW],  x3+DX3[NW],  NW);
      distributions->setDistributionInvForDirection(f[TE],  x1+DX1[BW],  x2+DX2[BW],  x3+DX3[BW],  BW);
      distributions->setDistributionInvForDirection(f[BE],  x1+DX1[TW],  x2+DX2[TW],  x3+DX3[TW],  TW);
      distributions->setDistributionInvForDirection(f[TNE], x1+DX1[BSW], x2+DX2[BSW], x3+DX3[BSW], BSW);
      distributions->setDistributionInvForDirection(f[TSE], x1+DX1[BNW], x2+DX2[BNW], x3+DX3[BNW], BNW);
      distributions->setDistributionInvForDirection(f[BNE], x1+DX1[TSW], x2+DX2[TSW], x3+DX3[TSW], TSW);
      distributions->setDistributionInvForDirection(f[BSE], x1+DX1[TNW], x2+DX2[TNW], x3+DX3[TNW], TNW);
      
	  h[E]   = htemp[E]   * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[E]   ;
	  h[NE]  = htemp[NE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[NE]  ;
	  h[SE]  = htemp[SE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[SE]  ;
	  h[TE]  = htemp[TE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[TE]  ;
	  h[BE]  = htemp[BE]  * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[BE]  ;
	  h[TNE] = htemp[TNE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[TNE] ;
	  h[TSE] = htemp[TSE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[TSE] ;
	  h[BNE] = htemp[BNE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[BNE] ;
	  h[BSE] = htemp[BSE] * (one_over_sqrt3 + vx1) + (1.0 - one_over_sqrt3 - vx1)*h[BSE] ;

	  distributionsH->setDistributionInvForDirection(h[E],   x1+DX1[W],   x2+DX2[W],   x3+DX3[W],   W);
	  distributionsH->setDistributionInvForDirection(h[NE],  x1+DX1[SW],  x2+DX2[SW],  x3+DX3[SW],  SW);
	  distributionsH->setDistributionInvForDirection(h[SE],  x1+DX1[NW],  x2+DX2[NW],  x3+DX3[NW],  NW);
	  distributionsH->setDistributionInvForDirection(h[TE],  x1+DX1[BW],  x2+DX2[BW],  x3+DX3[BW],  BW);
	  distributionsH->setDistributionInvForDirection(h[BE],  x1+DX1[TW],  x2+DX2[TW],  x3+DX3[TW],  TW);
	  distributionsH->setDistributionInvForDirection(h[TNE], x1+DX1[BSW], x2+DX2[BSW], x3+DX3[BSW], BSW);
	  distributionsH->setDistributionInvForDirection(h[TSE], x1+DX1[BNW], x2+DX2[BNW], x3+DX3[BNW], BNW);
	  distributionsH->setDistributionInvForDirection(h[BNE], x1+DX1[TSW], x2+DX2[TSW], x3+DX3[TSW], TSW);
	  distributionsH->setDistributionInvForDirection(h[BSE], x1+DX1[TNW], x2+DX2[TNW], x3+DX3[TNW], TNW);
	  
	  break;
   case W:
      f[W]   = ftemp[W]   * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[W]  ;
      f[NW]  = ftemp[NW]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[NW] ;
      f[SW]  = ftemp[SW]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[SW] ;
      f[TW]  = ftemp[TW]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[TW] ;
      f[BW]  = ftemp[BW]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[BW] ;
      f[TNW] = ftemp[TNW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[TNW];
      f[TSW] = ftemp[TSW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[TSW];
      f[BNW] = ftemp[BNW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[BNW];
      f[BSW] = ftemp[BSW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[BSW];

      distributions->setDistributionInvForDirection(f[W],   x1+DX1[E],   x2+DX2[E],   x3+DX3[E],     E);
      distributions->setDistributionInvForDirection(f[NW],  x1+DX1[SE],  x2+DX2[SE],  x3+DX3[SE],   SE);
      distributions->setDistributionInvForDirection(f[SW],  x1+DX1[NE],  x2+DX2[NE],  x3+DX3[NE],   NE);
      distributions->setDistributionInvForDirection(f[TW],  x1+DX1[BE],  x2+DX2[BE],  x3+DX3[BE],   BE);
      distributions->setDistributionInvForDirection(f[BW],  x1+DX1[TE],  x2+DX2[TE],  x3+DX3[TE],   TE);
      distributions->setDistributionInvForDirection(f[TNW], x1+DX1[BSE], x2+DX2[BSE], x3+DX3[BSE], BSE);
      distributions->setDistributionInvForDirection(f[TSW], x1+DX1[BNE], x2+DX2[BNE], x3+DX3[BNE], BNE);
      distributions->setDistributionInvForDirection(f[BNW], x1+DX1[TSE], x2+DX2[TSE], x3+DX3[TSE], TSE);
      distributions->setDistributionInvForDirection(f[BSW], x1+DX1[TNE], x2+DX2[TNE], x3+DX3[TNE], TNE);

	  h[W]   = htemp[W]  * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[W]  ;
	  h[NW]  = htemp[NW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[NW] ;
	  h[SW]  = htemp[SW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[SW] ;
	  h[TW]  = htemp[TW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[TW] ;
	  h[BW]  = htemp[BW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[BW] ;
	  h[TNW] = htemp[TNW]* (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[TNW];
	  h[TSW] = htemp[TSW]* (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[TSW];
	  h[BNW] = htemp[BNW]* (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[BNW];
	  h[BSW] = htemp[BSW]* (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*h[BSW];

	  distributionsH->setDistributionInvForDirection(h[W],   x1+DX1[E],   x2+DX2[E],   x3+DX3[E],     E);
	  distributionsH->setDistributionInvForDirection(h[NW],  x1+DX1[SE],  x2+DX2[SE],  x3+DX3[SE],   SE);
	  distributionsH->setDistributionInvForDirection(h[SW],  x1+DX1[NE],  x2+DX2[NE],  x3+DX3[NE],   NE);
	  distributionsH->setDistributionInvForDirection(h[TW],  x1+DX1[BE],  x2+DX2[BE],  x3+DX3[BE],   BE);
	  distributionsH->setDistributionInvForDirection(h[BW],  x1+DX1[TE],  x2+DX2[TE],  x3+DX3[TE],   TE);
	  distributionsH->setDistributionInvForDirection(h[TNW], x1+DX1[BSE], x2+DX2[BSE], x3+DX3[BSE], BSE);
	  distributionsH->setDistributionInvForDirection(h[TSW], x1+DX1[BNE], x2+DX2[BNE], x3+DX3[BNE], BNE);
	  distributionsH->setDistributionInvForDirection(h[BNW], x1+DX1[TSE], x2+DX2[TSE], x3+DX3[TSE], TSE);
	  distributionsH->setDistributionInvForDirection(h[BSW], x1+DX1[TNE], x2+DX2[TNE], x3+DX3[TNE], TNE);

      break;
   case N:
      f[N]   = ftemp[N]   * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[N]   ;
      f[NE]  = ftemp[NE]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[NE]  ;
      f[NW]  = ftemp[NW]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[NW]  ;
      f[TN]  = ftemp[TN]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[TN]  ;
      f[BN]  = ftemp[BN]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[BN]  ;
      f[TNE] = ftemp[TNE] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[TNE] ;
      f[TNW] = ftemp[TNW] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[TNW] ;
      f[BNE] = ftemp[BNE] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[BNE] ;
      f[BNW] = ftemp[BNW] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*f[BNW] ;

      distributions->setDistributionInvForDirection(f[N],   x1+DX1[S],   x2+DX2[S],   x3+DX3[S],     S);
      distributions->setDistributionInvForDirection(f[NE],  x1+DX1[SW],  x2+DX2[SW],  x3+DX3[SW],   SW);
      distributions->setDistributionInvForDirection(f[NW],  x1+DX1[SE],  x2+DX2[SE],  x3+DX3[SE],   SE);
      distributions->setDistributionInvForDirection(f[TN],  x1+DX1[BS],  x2+DX2[BS],  x3+DX3[BS],   BS);
      distributions->setDistributionInvForDirection(f[BN],  x1+DX1[TS],  x2+DX2[TS],  x3+DX3[TS],   TS);
      distributions->setDistributionInvForDirection(f[TNE], x1+DX1[BSW], x2+DX2[BSW], x3+DX3[BSW], BSW);
      distributions->setDistributionInvForDirection(f[TNW], x1+DX1[BSE], x2+DX2[BSE], x3+DX3[BSE], BSE);
      distributions->setDistributionInvForDirection(f[BNE], x1+DX1[TSW], x2+DX2[TSW], x3+DX3[TSW], TSW);
      distributions->setDistributionInvForDirection(f[BNW], x1+DX1[TSE], x2+DX2[TSE], x3+DX3[TSE], TSE);

	  h[N]   = htemp[N]   * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[N]   ;
	  h[NE]  = htemp[NE]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[NE]  ;
	  h[NW]  = htemp[NW]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[NW]  ;
	  h[TN]  = htemp[TN]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[TN]  ;
	  h[BN]  = htemp[BN]  * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[BN]  ;
	  h[TNE] = htemp[TNE] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[TNE] ;
	  h[TNW] = htemp[TNW] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[TNW] ;
	  h[BNE] = htemp[BNE] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[BNE] ;
	  h[BNW] = htemp[BNW] * (one_over_sqrt3 + vx2) + (1.0 - one_over_sqrt3 - vx2)*h[BNW] ;

	  distributionsH->setDistributionInvForDirection(h[N],   x1+DX1[S],   x2+DX2[S],   x3+DX3[S],     S);
	  distributionsH->setDistributionInvForDirection(h[NE],  x1+DX1[SW],  x2+DX2[SW],  x3+DX3[SW],   SW);
	  distributionsH->setDistributionInvForDirection(h[NW],  x1+DX1[SE],  x2+DX2[SE],  x3+DX3[SE],   SE);
	  distributionsH->setDistributionInvForDirection(h[TN],  x1+DX1[BS],  x2+DX2[BS],  x3+DX3[BS],   BS);
	  distributionsH->setDistributionInvForDirection(h[BN],  x1+DX1[TS],  x2+DX2[TS],  x3+DX3[TS],   TS);
	  distributionsH->setDistributionInvForDirection(h[TNE], x1+DX1[BSW], x2+DX2[BSW], x3+DX3[BSW], BSW);
	  distributionsH->setDistributionInvForDirection(h[TNW], x1+DX1[BSE], x2+DX2[BSE], x3+DX3[BSE], BSE);
	  distributionsH->setDistributionInvForDirection(h[BNE], x1+DX1[TSW], x2+DX2[TSW], x3+DX3[TSW], TSW);
	  distributionsH->setDistributionInvForDirection(h[BNW], x1+DX1[TSE], x2+DX2[TSE], x3+DX3[TSE], TSE);

      break;
   case S:
      f[S]   = ftemp[S]   * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[S]   ;
      f[SE]  = ftemp[SE]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[SE]  ;
      f[SW]  = ftemp[SW]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[SW]  ;
      f[TS]  = ftemp[TS]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[TS]  ;
      f[BS]  = ftemp[BS]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[BS]  ;
      f[TSE] = ftemp[TSE] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[TSE] ;
      f[TSW] = ftemp[TSW] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[TSW] ;
      f[BSE] = ftemp[BSE] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[BSE] ;
      f[BSW] = ftemp[BSW] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*f[BSW] ;

      distributions->setDistributionInvForDirection(f[S],   x1+DX1[N],   x2+DX2[N],   x3+DX3[N],     N);
      distributions->setDistributionInvForDirection(f[SE],  x1+DX1[NW],  x2+DX2[NW],  x3+DX3[NW],   NW);
      distributions->setDistributionInvForDirection(f[SW],  x1+DX1[NE],  x2+DX2[NE],  x3+DX3[NE],   NE);
      distributions->setDistributionInvForDirection(f[TS],  x1+DX1[BN],  x2+DX2[BN],  x3+DX3[BN],   BN);
      distributions->setDistributionInvForDirection(f[BS],  x1+DX1[TN],  x2+DX2[TN],  x3+DX3[TN],   TN);
      distributions->setDistributionInvForDirection(f[TSE], x1+DX1[BNW], x2+DX2[BNW], x3+DX3[BNW], BNW);
      distributions->setDistributionInvForDirection(f[TSW], x1+DX1[BNE], x2+DX2[BNE], x3+DX3[BNE], BNE);
      distributions->setDistributionInvForDirection(f[BSE], x1+DX1[TNW], x2+DX2[TNW], x3+DX3[TNW], TNW);
      distributions->setDistributionInvForDirection(f[BSW], x1+DX1[TNE], x2+DX2[TNE], x3+DX3[TNE], TNE);

	  h[S]   = htemp[S]   * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[S]   ;
	  h[SE]  = htemp[SE]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[SE]  ;
	  h[SW]  = htemp[SW]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[SW]  ;
	  h[TS]  = htemp[TS]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[TS]  ;
	  h[BS]  = htemp[BS]  * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[BS]  ;
	  h[TSE] = htemp[TSE] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[TSE] ;
	  h[TSW] = htemp[TSW] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[TSW] ;
	  h[BSE] = htemp[BSE] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[BSE] ;
	  h[BSW] = htemp[BSW] * (one_over_sqrt3 - vx2) + (1.0 - one_over_sqrt3 + vx2)*h[BSW] ;

	  distributionsH->setDistributionInvForDirection(h[S],   x1+DX1[N],   x2+DX2[N],   x3+DX3[N],     N);
	  distributionsH->setDistributionInvForDirection(h[SE],  x1+DX1[NW],  x2+DX2[NW],  x3+DX3[NW],   NW);
	  distributionsH->setDistributionInvForDirection(h[SW],  x1+DX1[NE],  x2+DX2[NE],  x3+DX3[NE],   NE);
	  distributionsH->setDistributionInvForDirection(h[TS],  x1+DX1[BN],  x2+DX2[BN],  x3+DX3[BN],   BN);
	  distributionsH->setDistributionInvForDirection(h[BS],  x1+DX1[TN],  x2+DX2[TN],  x3+DX3[TN],   TN);
	  distributionsH->setDistributionInvForDirection(h[TSE], x1+DX1[BNW], x2+DX2[BNW], x3+DX3[BNW], BNW);
	  distributionsH->setDistributionInvForDirection(h[TSW], x1+DX1[BNE], x2+DX2[BNE], x3+DX3[BNE], BNE);
	  distributionsH->setDistributionInvForDirection(h[BSE], x1+DX1[TNW], x2+DX2[TNW], x3+DX3[TNW], TNW);
	  distributionsH->setDistributionInvForDirection(h[BSW], x1+DX1[TNE], x2+DX2[TNE], x3+DX3[TNE], TNE);

      break;
   case T:
      f[T]   = ftemp[T]   * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[T]   ;
      f[TE]  = ftemp[TE]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[TE]  ;
      f[TW]  = ftemp[TW]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[TW]  ;
      f[TN]  = ftemp[TN]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[TN]  ;
      f[TS]  = ftemp[TS]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[TS]  ;
      f[TNE] = ftemp[TNE] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[TNE] ;
      f[TNW] = ftemp[TNW] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[TNW] ;
      f[TSE] = ftemp[TSE] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[TSE] ;
      f[TSW] = ftemp[TSW] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*f[TSW] ;

      distributions->setDistributionInvForDirection(f[T],   x1+DX1[B],   x2+DX2[B],   x3+DX3[B],     B);
      distributions->setDistributionInvForDirection(f[TE],  x1+DX1[BW],  x2+DX2[BW],  x3+DX3[BW],   BW);
      distributions->setDistributionInvForDirection(f[TW],  x1+DX1[BE],  x2+DX2[BE],  x3+DX3[BE],   BE);
      distributions->setDistributionInvForDirection(f[TN],  x1+DX1[BS],  x2+DX2[BS],  x3+DX3[BS],   BS);
      distributions->setDistributionInvForDirection(f[TS],  x1+DX1[BN],  x2+DX2[BN],  x3+DX3[BN],   BN);
      distributions->setDistributionInvForDirection(f[TNE], x1+DX1[BSW], x2+DX2[BSW], x3+DX3[BSW], BSW);
      distributions->setDistributionInvForDirection(f[TNW], x1+DX1[BSE], x2+DX2[BSE], x3+DX3[BSE], BSE);
      distributions->setDistributionInvForDirection(f[TSE], x1+DX1[BNW], x2+DX2[BNW], x3+DX3[BNW], BNW);
      distributions->setDistributionInvForDirection(f[TSW], x1+DX1[BNE], x2+DX2[BNE], x3+DX3[BNE], BNE);

	  h[T]   = htemp[T]   * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[T]   ;
	  h[TE]  = htemp[TE]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[TE]  ;
	  h[TW]  = htemp[TW]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[TW]  ;
	  h[TN]  = htemp[TN]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[TN]  ;
	  h[TS]  = htemp[TS]  * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[TS]  ;
	  h[TNE] = htemp[TNE] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[TNE] ;
	  h[TNW] = htemp[TNW] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[TNW] ;
	  h[TSE] = htemp[TSE] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[TSE] ;
	  h[TSW] = htemp[TSW] * (one_over_sqrt3 + vx3) + (1.0 - one_over_sqrt3 - vx3)*h[TSW] ;

	  distributionsH->setDistributionInvForDirection(h[T],   x1+DX1[B],   x2+DX2[B],   x3+DX3[B],     B);
	  distributionsH->setDistributionInvForDirection(h[TE],  x1+DX1[BW],  x2+DX2[BW],  x3+DX3[BW],   BW);
	  distributionsH->setDistributionInvForDirection(h[TW],  x1+DX1[BE],  x2+DX2[BE],  x3+DX3[BE],   BE);
	  distributionsH->setDistributionInvForDirection(h[TN],  x1+DX1[BS],  x2+DX2[BS],  x3+DX3[BS],   BS);
	  distributionsH->setDistributionInvForDirection(h[TS],  x1+DX1[BN],  x2+DX2[BN],  x3+DX3[BN],   BN);
	  distributionsH->setDistributionInvForDirection(h[TNE], x1+DX1[BSW], x2+DX2[BSW], x3+DX3[BSW], BSW);
	  distributionsH->setDistributionInvForDirection(h[TNW], x1+DX1[BSE], x2+DX2[BSE], x3+DX3[BSE], BSE);
	  distributionsH->setDistributionInvForDirection(h[TSE], x1+DX1[BNW], x2+DX2[BNW], x3+DX3[BNW], BNW);
	  distributionsH->setDistributionInvForDirection(h[TSW], x1+DX1[BNE], x2+DX2[BNE], x3+DX3[BNE], BNE);

      break;
   case B:
      f[B]   = ftemp[B]   * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[B]   ;
      f[BE]  = ftemp[BE]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[BE]  ;
      f[BW]  = ftemp[BW]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[BW]  ;
      f[BN]  = ftemp[BN]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[BN]  ;
      f[BS]  = ftemp[BS]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[BS]  ;
      f[BNE] = ftemp[BNE] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[BNE] ;
      f[BNW] = ftemp[BNW] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[BNW] ;
      f[BSE] = ftemp[BSE] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[BSE] ;
      f[BSW] = ftemp[BSW] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*f[BSW] ;

      distributions->setDistributionInvForDirection(f[B],   x1+DX1[T],   x2+DX2[T],   x3+DX3[T],     T);
      distributions->setDistributionInvForDirection(f[BE],  x1+DX1[TW],  x2+DX2[TW],  x3+DX3[TW],   TW);
      distributions->setDistributionInvForDirection(f[BW],  x1+DX1[TE],  x2+DX2[TE],  x3+DX3[TE],   TE);
      distributions->setDistributionInvForDirection(f[BN],  x1+DX1[TS],  x2+DX2[TS],  x3+DX3[TS],   TS);
      distributions->setDistributionInvForDirection(f[BS],  x1+DX1[TN],  x2+DX2[TN],  x3+DX3[TN],   TN);
      distributions->setDistributionInvForDirection(f[BNE], x1+DX1[TSW], x2+DX2[TSW], x3+DX3[TSW], TSW);
      distributions->setDistributionInvForDirection(f[BNW], x1+DX1[TSE], x2+DX2[TSE], x3+DX3[TSE], TSE);
      distributions->setDistributionInvForDirection(f[BSE], x1+DX1[TNW], x2+DX2[TNW], x3+DX3[TNW], TNW);
      distributions->setDistributionInvForDirection(f[BSW], x1+DX1[TNE], x2+DX2[TNE], x3+DX3[TNE], TNE);

	  h[B]   = htemp[B]   * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[B]   ;
	  h[BE]  = htemp[BE]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[BE]  ;
	  h[BW]  = htemp[BW]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[BW]  ;
	  h[BN]  = htemp[BN]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[BN]  ;
	  h[BS]  = htemp[BS]  * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[BS]  ;
	  h[BNE] = htemp[BNE] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[BNE] ;
	  h[BNW] = htemp[BNW] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[BNW] ;
	  h[BSE] = htemp[BSE] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[BSE] ;
	  h[BSW] = htemp[BSW] * (one_over_sqrt3 - vx3) + (1.0 - one_over_sqrt3 + vx3)*h[BSW] ;

	  distributionsH->setDistributionInvForDirection(h[B],   x1+DX1[T],   x2+DX2[T],   x3+DX3[T],     T);
	  distributionsH->setDistributionInvForDirection(h[BE],  x1+DX1[TW],  x2+DX2[TW],  x3+DX3[TW],   TW);
	  distributionsH->setDistributionInvForDirection(h[BW],  x1+DX1[TE],  x2+DX2[TE],  x3+DX3[TE],   TE);
	  distributionsH->setDistributionInvForDirection(h[BN],  x1+DX1[TS],  x2+DX2[TS],  x3+DX3[TS],   TS);
	  distributionsH->setDistributionInvForDirection(h[BS],  x1+DX1[TN],  x2+DX2[TN],  x3+DX3[TN],   TN);
	  distributionsH->setDistributionInvForDirection(h[BNE], x1+DX1[TSW], x2+DX2[TSW], x3+DX3[TSW], TSW);
	  distributionsH->setDistributionInvForDirection(h[BNW], x1+DX1[TSE], x2+DX2[TSE], x3+DX3[TSE], TSE);
	  distributionsH->setDistributionInvForDirection(h[BSE], x1+DX1[TNW], x2+DX2[TNW], x3+DX3[TNW], TNW);
	  distributionsH->setDistributionInvForDirection(h[BSW], x1+DX1[TNE], x2+DX2[TNE], x3+DX3[TNE], TNE);

      break;
   default:
      UB_THROW(UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
   }
}

