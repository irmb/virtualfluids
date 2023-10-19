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
//! \file SimpleSlipBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "SimpleSlipBCStrategy.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

SimpleSlipBCStrategy::SimpleSlipBCStrategy()
{
   BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
SimpleSlipBCStrategy::~SimpleSlipBCStrategy()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> SimpleSlipBCStrategy::clone()
{
   SPtr<BCStrategy> bc(new SimpleSlipBCStrategy());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void SimpleSlipBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void SimpleSlipBCStrategy::applyBC()
{
    using namespace vf::lbm::dir;

   real f[D3Q27System::ENDF+1];
   real feq[D3Q27System::ENDF+1];
   distributions->getDistributionInv(f, x1, x2, x3);
   real vx1, vx2, vx3, drho, rho;
   calcMacrosFct(f, drho, vx1, vx2, vx3);
   calcFeqFct(feq, drho, vx1, vx2, vx3);

   rho = vf::basics::constant::c1o1 + drho * compressibleFactor;

   UbTupleFloat3 normale = bcPtr->getNormalVector();
   real amp = vx1*val<1>(normale)+vx2*val<2>(normale)+vx3*val<3>(normale);

   vx1 = vx1 - amp * val<1>(normale); //normale zeigt von struktur weg!
   vx2 = vx2 - amp * val<2>(normale); //normale zeigt von struktur weg!
   vx3 = vx3 - amp * val<3>(normale); //normale zeigt von struktur weg!

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fdir];
         real velocity = vf::basics::constant::c0o1;
         switch (invDir)
         {
         case DIR_P00: velocity = (vf::basics::constant::c4o9*(+vx1)); break;      //(2/cs^2)(=6)*rho_0(=1 bei imkompr)*wi*u*ei mit cs=1/sqrt(3)
         case DIR_M00: velocity = (vf::basics::constant::c4o9*(-vx1)); break;      //z.B. aus paper manfred MRT LB models in three dimensions (2002)   
         case DIR_0P0: velocity = (vf::basics::constant::c4o9*(+vx2)); break;
         case DIR_0M0: velocity = (vf::basics::constant::c4o9*(-vx2)); break;
         case DIR_00P: velocity = (vf::basics::constant::c4o9*(+vx3)); break;
         case DIR_00M: velocity = (vf::basics::constant::c4o9*(-vx3)); break;
         case DIR_PP0: velocity = (vf::basics::constant::c1o9*(+vx1+vx2)); break;
         case DIR_MM0: velocity = (vf::basics::constant::c1o9*(-vx1-vx2)); break;
         case DIR_PM0: velocity = (vf::basics::constant::c1o9*(+vx1-vx2)); break;
         case DIR_MP0: velocity = (vf::basics::constant::c1o9*(-vx1+vx2)); break;
         case DIR_P0P: velocity = (vf::basics::constant::c1o9*(+vx1+vx3)); break;
         case DIR_M0M: velocity = (vf::basics::constant::c1o9*(-vx1-vx3)); break;
         case DIR_P0M: velocity = (vf::basics::constant::c1o9*(+vx1-vx3)); break;
         case DIR_M0P: velocity = (vf::basics::constant::c1o9*(-vx1+vx3)); break;
         case DIR_0PP: velocity = (vf::basics::constant::c1o9*(+vx2+vx3)); break;
         case DIR_0MM: velocity = (vf::basics::constant::c1o9*(-vx2-vx3)); break;
         case DIR_0PM: velocity = (vf::basics::constant::c1o9*(+vx2-vx3)); break;
         case DIR_0MP: velocity = (vf::basics::constant::c1o9*(-vx2+vx3)); break;
         case DIR_PPP: velocity = (vf::basics::constant::c1o36*(+vx1+vx2+vx3)); break;
         case DIR_MMM: velocity = (vf::basics::constant::c1o36*(-vx1-vx2-vx3)); break;
         case DIR_PPM: velocity = (vf::basics::constant::c1o36*(+vx1+vx2-vx3)); break;
         case DIR_MMP: velocity = (vf::basics::constant::c1o36*(-vx1-vx2+vx3)); break;
         case DIR_PMP: velocity = (vf::basics::constant::c1o36*(+vx1-vx2+vx3)); break;
         case DIR_MPM: velocity = (vf::basics::constant::c1o36*(-vx1+vx2-vx3)); break;
         case DIR_PMM: velocity = (vf::basics::constant::c1o36*(+vx1-vx2-vx3)); break;
         case DIR_MPP: velocity = (vf::basics::constant::c1o36*(-vx1+vx2+vx3)); break;
         default: throw UbException(UB_EXARGS, "unknown error");
         }
         real fReturn = f[invDir] - velocity * rho;
         distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }
}