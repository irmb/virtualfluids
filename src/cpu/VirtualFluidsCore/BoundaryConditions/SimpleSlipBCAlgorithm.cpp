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
//! \file SimpleSlipBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "SimpleSlipBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

SimpleSlipBCAlgorithm::SimpleSlipBCAlgorithm()
{
   BCAlgorithm::type = BCAlgorithm::SimpleSlipBCAlgorithm;
   BCAlgorithm::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
SimpleSlipBCAlgorithm::~SimpleSlipBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> SimpleSlipBCAlgorithm::clone()
{
   SPtr<BCAlgorithm> bc(new SimpleSlipBCAlgorithm());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void SimpleSlipBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void SimpleSlipBCAlgorithm::applyBC()
{
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal feq[D3Q27System::ENDF+1];
   distributions->getDistributionInv(f, x1, x2, x3);
   LBMReal vx1, vx2, vx3, drho, rho;
   calcMacrosFct(f, drho, vx1, vx2, vx3);
   calcFeqFct(feq, drho, vx1, vx2, vx3);

   rho = 1.0 + drho * compressibleFactor;

   UbTupleFloat3 normale = bcPtr->getNormalVector();
   LBMReal amp = vx1*val<1>(normale)+vx2*val<2>(normale)+vx3*val<3>(normale);

   vx1 = vx1 - amp * val<1>(normale); //normale zeigt von struktur weg!
   vx2 = vx2 - amp * val<2>(normale); //normale zeigt von struktur weg!
   vx3 = vx3 - amp * val<3>(normale); //normale zeigt von struktur weg!

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fdir];
         LBMReal velocity = 0.0;
         switch (invDir)
         {
         case D3Q27System::DIR_P00: velocity = (UbMath::c4o9*(+vx1)); break;      //(2/cs^2)(=6)*rho_0(=1 bei imkompr)*wi*u*ei mit cs=1/sqrt(3)
         case D3Q27System::DIR_M00: velocity = (UbMath::c4o9*(-vx1)); break;      //z.B. aus paper manfred MRT LB models in three dimensions (2002)   
         case D3Q27System::DIR_0P0: velocity = (UbMath::c4o9*(+vx2)); break;
         case D3Q27System::DIR_0M0: velocity = (UbMath::c4o9*(-vx2)); break;
         case D3Q27System::DIR_00P: velocity = (UbMath::c4o9*(+vx3)); break;
         case D3Q27System::DIR_00M: velocity = (UbMath::c4o9*(-vx3)); break;
         case D3Q27System::DIR_PP0: velocity = (UbMath::c1o9*(+vx1+vx2)); break;
         case D3Q27System::DIR_MM0: velocity = (UbMath::c1o9*(-vx1-vx2)); break;
         case D3Q27System::DIR_PM0: velocity = (UbMath::c1o9*(+vx1-vx2)); break;
         case D3Q27System::DIR_MP0: velocity = (UbMath::c1o9*(-vx1+vx2)); break;
         case D3Q27System::DIR_P0P: velocity = (UbMath::c1o9*(+vx1+vx3)); break;
         case D3Q27System::DIR_M0M: velocity = (UbMath::c1o9*(-vx1-vx3)); break;
         case D3Q27System::DIR_P0M: velocity = (UbMath::c1o9*(+vx1-vx3)); break;
         case D3Q27System::DIR_M0P: velocity = (UbMath::c1o9*(-vx1+vx3)); break;
         case D3Q27System::DIR_0PP: velocity = (UbMath::c1o9*(+vx2+vx3)); break;
         case D3Q27System::DIR_0MM: velocity = (UbMath::c1o9*(-vx2-vx3)); break;
         case D3Q27System::DIR_0PM: velocity = (UbMath::c1o9*(+vx2-vx3)); break;
         case D3Q27System::DIR_0MP: velocity = (UbMath::c1o9*(-vx2+vx3)); break;
         case D3Q27System::DIR_PPP: velocity = (UbMath::c1o36*(+vx1+vx2+vx3)); break;
         case D3Q27System::DIR_MMM: velocity = (UbMath::c1o36*(-vx1-vx2-vx3)); break;
         case D3Q27System::DIR_PPM: velocity = (UbMath::c1o36*(+vx1+vx2-vx3)); break;
         case D3Q27System::DIR_MMP: velocity = (UbMath::c1o36*(-vx1-vx2+vx3)); break;
         case D3Q27System::DIR_PMP: velocity = (UbMath::c1o36*(+vx1-vx2+vx3)); break;
         case D3Q27System::DIR_MPM: velocity = (UbMath::c1o36*(-vx1+vx2-vx3)); break;
         case D3Q27System::DIR_PMM: velocity = (UbMath::c1o36*(+vx1-vx2-vx3)); break;
         case D3Q27System::DIR_MPP: velocity = (UbMath::c1o36*(-vx1+vx2+vx3)); break;
         default: throw UbException(UB_EXARGS, "unknown error");
         }
         LBMReal fReturn = f[invDir] - velocity * rho;
         distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }
}